"""Fast Structured Grid Navier-Stokes solver.
"""
import numpy
import os
import FastC.fastc
from . import fasts
from . import FastS
__version__ = FastS.__version__


try:
    import Converter.PyTree as C
    import Converter.Mpi as Cmpi
    import Converter.Internal as Internal
    import Connector
    import Connector.PyTree as X
    import Connector.OversetData as XOD
    import FastC.PyTree as FastC
    import math
    import time as Time
    import Connector.Mpi as Xmpi
    import RigidMotion.PyTree as R

    #import timeit
    #import KCore.Dist as Dist
except:
  raise ImportError("FastS: requires Converter, Connector, Fast modules.")

try:
    OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS']
    OMP_NUM_THREADS = int(OMP_NUM_THREADS)
except: OMP_NUM_THREADS = 1

try: range = xrange
except: pass

# Variable alignement pour vectorisation
#CACHELINE = Dist.getCacheLine()

#==============================================================================
# generation maillage pour tble
#==============================================================================
def _stretch(coord, nbp, x1, x2, dx1, dx2, ityp):
    fasts._stretch(coord, nbp, x1, x2, dx1, dx2, ityp)
    return None

#==============================================================================
# compute in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None, tc2=None, graph2=None, layer="c", NIT=1, ucData=None, vtune=None):
    """Compute a given number of iterations."""
    gradP      =False
    TBLE       =False
    isWireModel=False
    if tc is not None:
        base       = Internal.getBases(tc)[0]
        solverIBC  = Internal.getNodeFromName(base ,'.Solver#IBCdefine')
        if solverIBC is not None:
            gradP      = eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isgradP')))
            TBLE       = eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isTBLE')))
            isWireModel= eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isWireModel')))

    if isWireModel and tc2 is not None:
        print("Wire model does not currently support tc2 options...exiting...")
        exit()
        
    bases  = Internal.getNodesFromType1(t , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud

    zones = Internal.getZones(t)

    ## Note: This error is currently placed as a peculiar behavior was observed for a specific
    ##       test case.
    ##       For identical parameters, mesh, running conditions, etc. the "c layer" mode resulted
    ##       in wrong results while the "python layer" gave the correct results.
    ##       However, for simplified regression & verification test cases "c layer" and "python layer"
    ##       resulted in the same results. Therefore the below warning is just for precaution and if needed can be
    ##       commented out.
        
    d = Internal.getNodeFromName1(t, '.Solver#define')
    a = Internal.getNodeFromName1(d, 'temporal_scheme')
    if Internal.getValue(a) == 'explicit_local' and layer=="c":
        raise ValueError("compute: 'explicit_local' & 'c layer' is not implemented.")

    node = Internal.getNodeFromName1(t, '.Solver#define')
    omp_node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode  = FastC.OMP_MODE
    if omp_node is not None: ompmode = Internal.getValue(omp_node)

    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])

    #### a blinder...
    param_int_firstZone = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1]
    itypcp = param_int_firstZone[29]
    rk     = param_int_firstZone[52]
    exploc = param_int_firstZone[54]
    #### a blinder...
    if nitrun == 1: print('Info: using layer trans=%s (ompmode=%d)'%(layer, ompmode))

    tps_cp=  Time.time(); tps_cp=  tps_cp-tps_cp     
    tps_tr=  Time.time(); tps_tr=  tps_tr-tps_tr  
    if layer == "Python":
        
      if exploc==1 and tc is not None:        
         tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
         if tc_compact is not None:
                param_real_tc= tc_compact[1]
                param_int_tc = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]

                zones_tc    = Internal.getZones(tc)
                nbcomIBC    = param_int_tc[1]
                shift_graph = nbcomIBC + param_int_tc[2+nbcomIBC] + 2
                comm_P2P    = param_int_tc[0]
                pt_ech      = param_int_tc[comm_P2P + shift_graph]
                dest        = param_int_tc[pt_ech]

      hookTransfer = []        

      for nstep in range(1, nitmax+1): # pas RK ou ssiterations
         hook1 = FastC.HOOK.copy()
         hook1.update(  FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, 0) )

         skip = 0
         if hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1: skip = 1

         # calcul Navier_stokes + appli CL
         if skip == 0:
            nstep_deb = nstep
            nstep_fin = nstep
            layer_mode= 0
            nit_c     = 1
            #t0=Time.time()
            tic=Time.time()
            fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c, hook1)
            tps_cp +=Time.time()-tic 
            #print('t_compute = %f'%(Time.time() - t0))

            timelevel_target = int(dtloc[4])

            # dtloc GJeanmasson
            if exploc==1 and tc is not None:
               fasts.dtlocal2para_(zones, zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1, dest)
               
               if    nstep%2 == 0 and itypcp == 2: vars = ['Density'  ]
               elif  nstep%2 == 1 and itypcp == 2: vars = ['Density_P1']
               _applyBC(zones,metrics, hook1, nstep, var=vars[0])
               
               FastC.switchPointers2__(zones,nitmax,nstep)

               # Ghostcell
               if    nstep%2 == 0 and itypcp == 2: vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
               elif  nstep%2 == 1 and itypcp == 2: vars = ['Density_P1']
               _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, nitmax=nitmax, rk=rk, exploc=exploc,isWireModel=isWireModel)

               fasts.recup3para_(zones,zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1)

               if nstep%2 == 0:
                   vars = ['Density']
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep,  ompmode, hook1, nitmax=nitmax, rk=rk, exploc=exploc, num_passage= 2,isWireModel=isWireModel)

               if   nstep%2 == 0 and itypcp == 2: vars = ['Density'  ]
               elif nstep%2 == 1 and itypcp == 2: vars = ['Density_P1']
               _applyBC(zones, metrics, hook1, nstep, var=vars[0])
               
            else:
              #Ghostcell
              vars = FastC.varsP
              if nstep%2 == 0 and itypcp == 2: vars = FastC.varsN  # Choix du tableau pour application transfer et BC
              #t0=Time.time()
              tic=Time.time()

              if not tc2:
                _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, TBLE=TBLE, gradP=gradP,isWireModel=isWireModel)
              else:
                _fillGhostcells2(zones, tc, tc2, metrics, timelevel_target, vars, nstep, ompmode, hook1, TBLE=TBLE, gradP=gradP)
                
              #print('t_fillGhost = ',  Time.time() - t0 ,'nstep =', nstep)

              # Add unsteady Chimera transfers here
              if ucData is not None:
                (graphX, intersectionDict, dictOfADT,
                dictOfNobOfRcvZones, dictOfNozOfRcvZones,
                dictOfNobOfDnrZones, dictOfNozOfDnrZones,
                dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC,
                time, procDict, interpInDnrFrame, varType, tfreq, order, verbose) = ucData
                
                if nstep%2 == 0 and itypcp == 2: 
                    VARS = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
                    if varType == 1: VARS += ['TurbulentSANuTilde']
                else: 
                    VARS = ['Density_P1', 'VelocityX_P1', 'VelocityY_P1', 'VelocityZ_P1', 'Temperature_P1']
                    if varType == 1: VARS += ['TurbulentSANuTilde_P1']
                for v in VARS: C._cpVars(t, 'centers:'+v, tc, v)
                C._cpVars(t, "centers:cellN", tc, "cellN")

                if nstep == nitmax or nstep%tfreq == 0:
                    Xmpi._transfer2(t, tc, VARS, graphX, intersectionDict, dictOfADT, 
                                    dictOfNobOfRcvZones, dictOfNozOfRcvZones,
                                    dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
                                    dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC, 
                                    time=time, absFrame=True,
                                    procDict=procDict, cellNName='cellN#Motion', 
                                    interpInDnrFrame=interpInDnrFrame, order=order, 
                                    hook=hookTransfer, verbose=verbose)
                    #print('t_transfert = ',  Time.time() - t0 ,'nstep =', nstep)
              #print('t_transferts = %f'%(Time.time() - t0)
              tps_tr += Time.time()-tic

    else: ### layer C
      if isWireModel:
          raise ValueError("compute: FastS.PyTree.py C layer doesn't currently support the wire mesh model.")

      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT

      if tc is None:
          FastC.HOOK['param_int_tc']  = None
          FastC.HOOK['param_real_tc'] = None
      tic=Time.time()
      tps_tr =fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c, FastC.HOOK)
      tps_cp +=Time.time()-tic - tps_tr 
    # switch pointer a la fin du pas de temps
    if exploc==1 and tc is not None:
         if layer == 'Python': FastC.switchPointers__(zones, 1, 3)
         else: FastC.switchPointers3__(zones, nitmax)
    else:
         case = NIT%3
         if case != 0 and itypcp < 2: FastC.switchPointers__(zones, case)
         if case != 0 and itypcp ==2: FastC.switchPointers__(zones, case, nitmax%2+2)
    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT          = 1

    if vtune is None:
      return None
    else:
      return tps_cp,tps_tr

#==============================================================================
# compute in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute_matvec(t, metrics, no_vect_test, tc=None, graph=None):

    own   = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud

    zones = Internal.getZones(t)

    node = Internal.getNodeFromName1(t, '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = FastC.OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    nstep = 1

    hook1       = FastC.HOOK.copy()
    nitrun      = 0
    hook1.update(  FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, 0) )

    fasts._matvecPT(zones, metrics, nitrun, no_vect_test, ompmode, hook1)

    return None

#==============================================================================
# alloue retourne la metrique
#==============================================================================
def allocate_metric(t):

    own          = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc        = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
    dtloc_numpy  = Internal.getValue(dtloc)
    nssiter      = int(dtloc_numpy[0])

    zones  = Internal.getZones(t)
    metrics=[]
    for z in zones:
        motion ='none'
        b = Internal.getNodeFromName2(z, 'motion')
        if b is not None: motion = Internal.getValue(b)
        num = Internal.getNodeFromName1(z, '.Solver#ownData')
        if num is None:
            raise ValueError("metric: numerics is missing for zone %s."%z[0])
        if motion == 'rigid' or motion == 'deformation':
            grids = Internal.getNodesFromType1(z, 'GridCoordinates_t')
            if len(grids) == 1:
               grid_init = Internal.copyTree(grids[0])
               grid_init[0] = 'GridCoordinates#Init'
               Internal.addChild(z, grid_init, pos=-1) # last
        metrics.append(fasts.allocate_metric(z, nssiter))
    return metrics

#==============================================================================
# alloue retourne tableau ssor
#==============================================================================
def allocate_ssor(t, metrics, hook, ompmode):
    own         = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc       = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
    dtloc_numpy = Internal.getValue(dtloc)
    nssiter = int(dtloc_numpy[0])

    zones = Internal.getZones(t)
    ssors = []
    ssors = fasts.allocate_ssor(zones, metrics, nssiter, hook, ompmode)

    return ssors

#==============================================================================
# alloue retourne la metrique
#==============================================================================
def _init_metric(t, metrics, hook):
    zones        = Internal.getZones(t)

    FastC.fastc.init_metric(zones, metrics, hook)

    c=0
    for metric in metrics:
       z = zones[c]
       param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
       if param_int[34]!= 2: del metric[7]  # on efface l'info sur maille degen si pas besoin de recalculer metric mesh deforme
       c+=1

    return None

#==============================================================================
# Initialisation parametre calcul: calcul metric + var primitive + compactage
# + alignement + placement DRAM
#==============================================================================
def warmup(t, tc, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None, verbose=0):
    """Perform necessary operations for the solver to run."""
    Re  =-1
    Lref= 1.
    gradP      =False
    isWireModel=False
    if tc is not None:
        base       = Internal.getBases(tc)[0]
        solverIBC  = Internal.getNodeFromName(base ,'.Solver#IBCdefine')
        if solverIBC is not None:
            gradP      = eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isgradP')))
            isWireModel= eval(Internal.getValue(Internal.getNodeFromName(solverIBC, 'isWireModel')))

            Re  = Internal.getValue(Internal.getNodeFromName(solverIBC, 'Reref'))
            Lref= Internal.getValue(Internal.getNodeFromName(solverIBC, 'Lref'))
            
    # compute info linelets
    nbpts_linelets = 0
    if tc is not None:
      if Re > 0: #Adaptive method
        h0, hn, nbpts_linelets = computeLineletsInfo(tc, Re=Re, Lref=Lref, q=1.1)
      elif Re == 0: 
        h0, hn, nbpts_linelets = computeLineletsInfo2(tc, q=1.1)
      else: #Alferez' og method
        h0, nbpts_linelets = 1.e-6, 45

    first = Internal.getNodeFromName1(t, 'NbptsLinelets')
    if first is None: Internal.createUniqueChild(t, 'NbptsLinelets', 'DataArray_t', value=nbpts_linelets)
    
    # Get omp_mode
    ompmode = FastC.OMP_MODE
    node = Internal.getNodeFromName1(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastC._reorder(t, tc, ompmode)
    # Construction param_int et param_real des zones
    FastC._buildOwnData(t, Padding)
    t = Internal.rmNodesByName(t, 'NbptsLinelets')
     
    #init hook necessaire pour info omp
    tmp     = Internal.getNodeFromName1(t, '.Solver#ownData')
    dtlocPy = Internal.getNodeFromName1(tmp, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtlocPy)                       # tab numpy
    zones = Internal.getZones(t)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it
    # allocation d'espace dans param_int pour stockage info openmp  
    FastC._build_omp(t)
    # alloue metric: tijk, ventijk, ssiter_loc
    # init         : ssiter_loc
    metrics = allocate_metric(t)
    # Contruction BC_int et BC_real pour CL
    FastC._BCcompact(t)
    # Contruction Conservative_int
    FastC._Fluxcompact(t) 
    #determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)
    for nstep in range(1, int(dtloc[0])+1):
        hook1       = FastC.HOOK.copy()
        hook1.update(  FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, 1, nstep, verbose) )   

    #init metric
    _init_metric(t, metrics, hook1)

    ssors = allocate_ssor(t, metrics, hook1, ompmode)

    # compact + align + init numa
    rmConsVars = True
    adjoint = Adjoint

    t, FastC.FIRST_IT, zones2compact = FastC.createPrimVars(t, ompmode, rmConsVars, adjoint, gradP,isWireModel)
    FastC.HOOK['FIRST_IT'] = FastC.FIRST_IT

    zones = Internal.getZones(t)

    #compactage des champs en fonction option de calcul
    count = -1
    if ompmode == 1: count = 0
    for data in zones2compact:
        if ompmode == 1: count += 1
        zone    = data[0]
        varnames= data[1]
        for fields in varnames:
            FastC._compact(zone, fields=fields, mode=count, dtloc=dtlocPy)

    #corection pointeur ventijk si ale=0: pointeur Ro perdu par compact.
    c   = 0
    ale = 0
    for z in zones:
        motion = 'none'
        b = Internal.getNodeFromName2(z, 'motion')
        if b is not None: motion = Internal.getValue(b)
        if motion == 'none':
            sol = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
            ro = Internal.getNodeFromName1(sol, 'Density')
            metrics[c][2] = ro[1]
        elif motion == 'deformation': ale = 2
        else: ale = 1
        c += 1


    # mise a jour vitesse entrainememnt
    #
    #t0=timeit.default_timer()
    if ale == 1 and infos_ale is not None:
        print("ale actif. Teta et tetap=", infos_ale)
        teta = infos_ale[0]; tetap = infos_ale[1]
        FastC._motionlaw(t, teta, tetap)
        _computeVelocityAle(t, metrics)
    elif ale == 2:
        first = Internal.getNodeFromName1(t, 'Time')
        if first is not None: time = Internal.getValue(first)
        else: time = 0.
        R._evalPosition(t, time)
        R._evalGridSpeed(t, time)
        copy_velocity_ale(t, metrics)

    #t1=timeit.default_timer()
    #print("cout mise a jour vitesse entr= ", t1-t0)
    #
    # Compactage arbre transfert et mise a jour FastC.HOOK
    #

    if tc is not None:
       # Add linelets arrays in the FastC.HOOK for ODE-based WallModel (IBC)
       # nbpts_linelets = 0
       # _createTBLESA(tc,nbpts_linelets)

       if Re > 0: #Adaptive method
         # h0, hn, nbpts_linelets = computeLineletsInfo(tc, Re=Re, Lref=Lref, q=1.1)
         _createTBLESA(tc, h0=h0, hn=hn, nbpts_linelets=nbpts_linelets)
         _createTBLESA2(t, tc, h0=h0, hn=hn, nbpts_linelets=nbpts_linelets)
       elif Re == 0: 
         # h0, hn, nbpts_linelets = computeLineletsInfo2(tc, q=1.1)
         _createTBLESA(tc, h0=h0, hn=hn, nbpts_linelets=nbpts_linelets)
         _createTBLESA2(t, tc, h0=h0, hn=hn, nbpts_linelets=nbpts_linelets)
       else: #Alferez' og method
         # h0, nbpts_linelets = 1.e-6, 45
         _createTBLESA(tc, h0=h0, hn=-1, nbpts_linelets=nbpts_linelets)
         _createTBLESA2(t, tc, h0=h0, hn=-1, nbpts_linelets=nbpts_linelets)
 
       X.miseAPlatDonorTree__(zones, tc, graph=graph, list_graph=list_graph, nbpts_linelets=nbpts_linelets)

       FastC.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
       param_real_tc = Internal.getNodeFromName1(tc, 'Parameter_real')
       if param_real_tc is not None: FastC.HOOK['param_real_tc']= param_real_tc[1]

       # FastC.HOOK['param_real_tc'][ 58 ] = nbpts_linelets

       # # Add linelets arrays in the FastC.HOOK for ODE-based WallModel (IBC)
       # nbpts_linelets = 45
       # _createTBLESA(tc,nbpts_linelets)
       # _createTBLESA2(tc,nbpts_linelets)

    else:
        FastC.HOOK['param_real_tc'] = None
        FastC.HOOK['param_int_tc']  = None

    if ssors is not []:
        FastC.HOOK['ssors'] = ssors
    else:
        FastC.HOOK['ssors'] = None

    # Compactage arbre moyennes stat
    #
    if tmy is not None:
        sol = Internal.getNodesFromName3(tmy, 'FlowSolution#Centers')
        var = Internal.getNodesFromType1(sol[0], 'DataArray_t')
        varmy = []
        for v in var: varmy.append('centers:'+v[0])
        FastC._compact(tmy, fields=varmy)

    #
    # remplissage ghostcells
    #
    hook1["lexit_lu"]= 0
    nstep            = 1
    nitrun           = 0

    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4])
    
    _fillGhostcells(zones, tc, metrics, timelevel_target, ['Density'], nstep,  ompmode, hook1,isWireModel=isWireModel)

    if tc is not None: C._rmVars(tc, 'FlowSolution')
    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    fasts._computePT_mut(zones, metrics, hook1)

    return (t, tc, metrics)

#==============================================================================
# Prepare IBC ODE (TBLE + Spalart 1D)
# Info stored and compacted in linelets_int and linelets_real
# Seq only
#==============================================================================
def _createTBLESA(tc, h0, hn, nbpts_linelets=45):
      nfields_lin  = 6 # u,nutild,y,matm,mat,matp
      count_racIBC = 0
      addr_IBC = 0
      addrIBC  = []
      size_IBC = 0

      linelets_int_tc  = Internal.getNodeFromName(tc,'linelets_int')
      linelets_real_tc = Internal.getNodeFromName(tc,'linelets_real')
      if linelets_int_tc is not None:
         FastC.HOOK['linelets_int'] = Internal.getValue(linelets_int_tc)
      if linelets_real_tc is not None:
         FastC.HOOK['linelets_real']= Internal.getValue(linelets_real_tc)

      zones_tc = Internal.getZones(tc)

      for z in zones_tc:
          subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
          if subRegions is not None:
                for s in subRegions:
                  zsrname = s[0]
                  sname  = zsrname[0:2]
                  if sname == 'IB':
                     zsrname = zsrname.split('_')
                     if len(zsrname)<3:
                        print('Warning: createTBLESA: non consistent with the version of IBM preprocessing.')
                     else:
                       if zsrname[1] == '6':
                        print('Using TBLE SA wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        # size_IBC      = size_IBC + Nbpts_D*nbpts_linelets*nfields_lin
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                       elif zsrname[1] == '15':
                        print('Using MuskerSA wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        # size_IBC      = size_IBC + Nbpts_D*nbpts_linelets*nfields_lin
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1

      # addrIBC.append(size_IBC)
      # print "Nb Rac =",count_racIBC-1

      if size_IBC > 0:
        if FastC.HOOK['linelets_int'] is None:   # No TBLE Data in tc
         print('TBLE data missing in tc, they will be built...')

         # Nbpts_Dtot = size_IBC/(nbpts_linelets*nfields_lin)
         Nbpts_Dtot = size_IBC//(nbpts_linelets*nfields_lin+1)
         _addrIBC = numpy.asarray([nbpts_linelets,count_racIBC-1,Nbpts_Dtot] + addrIBC + [0]*Nbpts_Dtot ,dtype=numpy.int32)

         FastC.HOOK['linelets_int'] = _addrIBC

         linelets = numpy.zeros(size_IBC + Nbpts_Dtot, dtype=numpy.float64)


         # ====================== Debug (Interp linelets from fine grid) ======================

         # velxData = Internal.getNodeFromName(t,'VelocityX')
         # velxDataval = Internal.getValue(velxData)

         # turbSAData    = Internal.getNodeFromName(t,'TurbulentSANuTilde')
         # turbSADataval = Internal.getValue(turbSAData)


         # coordYdata = Internal.getNodeFromName(t,'CoordinateY')
         # coordYdataval = Internal.getValue(coordYdata)
         # Iw = numpy.argwhere(coordYdataval<0)
         # indexWall = Iw[-1][1]


         # ptlist = Internal.getNodeFromName(tc,'PointListDonor')#'TurbulentSANuTilde')#'VelocityX')
         # ptlistval = Internal.getValue(ptlist)

         # dimI   = Internal.getZoneDim(zones[0])[1]

         # nbptD = 190
         # ycoord = numpy.zeros(nbptD,dtype=numpy.float64)
         # field  = numpy.zeros(nbptD,dtype=numpy.float64)


         shift = 0
         for z in zones_tc:
                       subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
                       if subRegions is not None:
                          for s in subRegions:
                              zsrname = s[0]
                              sname  = zsrname[0:2]
                              if sname == 'IB':
                                 zsrname = zsrname.split('_')
                                 if len(zsrname)<3:
                                   print('Warning: createTBLESA: non consistent with the version of IBM preprocessing.')
                                 else:
                                   if zsrname[1] == '6' or zsrname[1] == '15':
                                     pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                                     Nbpts_D       = numpy.shape(pointlistD[1])[0]

                                     zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')
                                     valzI = Internal.getValue(zI)
                                     yI = Internal.getNodeFromName1(s,'CoordinateY_PI')
                                     valyI = Internal.getValue(yI)
                                     xI = Internal.getNodeFromName1(s,'CoordinateX_PI')
                                     valxI = Internal.getValue(xI)
                                     zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')
                                     valzW = Internal.getValue(zW)
                                     yW = Internal.getNodeFromName1(s,'CoordinateY_PW')
                                     valyW = Internal.getValue(yW)
                                     xW = Internal.getNodeFromName1(s,'CoordinateX_PW')
                                     valxW = Internal.getValue(xW)

                                     if Nbpts_D > 1:

                                         for noind in range(0, Nbpts_D):

                                             b0 = valxI[noind]-valxW[noind]
                                             b1 = valyI[noind]-valyW[noind]
                                             b2 = valzI[noind]-valzW[noind]
                                             normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)

                                             if hn < 0:
                                               hn = 0.6*normb

                                             # _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                             _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets], nbpts_linelets, 0.0, normb, h0, hn, 2)
                                     else:

                                             b0 = valxI-valxW
                                             b1 = valyI-valyW
                                             b2 = valzI-valzW
                                             normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)

                                             if hn < 0:
                                               hn = 0.6*normb

                                             # _stretch(linelets[shift : shift + nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                             _stretch(linelets[shift : shift + nbpts_linelets], nbpts_linelets, 0.0, normb, h0, hn, 2)


                                     #Internal.createUniqueChild(s, 'CoordinateN_ODE', 'DataArray_t',
                                     #                           linelets[shift : shift + Nbpts_D*nbpts_linelets])
                                     #Internal.createUniqueChild(s, 'VelocityT_ODE', 'DataArray_t',
                                     #                           linelets[shift + Nbpts_D*nbpts_linelets : shift + 2*Nbpts_D*nbpts_linelets ])
                                     #Internal.createUniqueChild(s, 'TurbulentSANuTilde_ODE', 'DataArray_t',
                                     #                           linelets[shift + 2*Nbpts_D*nbpts_linelets  : shift + 3*Nbpts_D*nbpts_linelets ])


                                     # ========================= Debug interp from fine grid =============

                                     # for noind in range(0,Nbpts_D):
                                     #     ltarget = ptlistval[noind]
                                     #     itarget = ltarget%(dimI-1)

                                     #     ycoord[:] = coordYdataval[itarget,indexWall+1:indexWall+nbptD+1,0]
                                     #     field[:]  =   velxDataval[itarget,indexWall+1:indexWall+nbptD+1,0] # VelocityX

                                     #     fasts._interpfromzone(nbptD,nbpts_linelets,ycoord,
                                     #                                       linelets[shift : shift + nbpts_linelets],
                                     #                                       field,
                                     #                                       linelets[shift + Nbpts_D*nbpts_linelets  + noind*nbpts_linelets : shift + Nbpts_D*nbpts_linelets + (noind+1)*nbpts_linelets])

                                     #     field[:]  =   turbSADataval[itarget,indexWall+1:indexWall+nbptD+1,0] # TurbulentSANutTilde

                                     #     fasts._interpfromzone(nbptD,nbpts_linelets,ycoord,
                                     #                                       linelets[shift : shift + nbpts_linelets],
                                     #                                       field,
                                     #                                       linelets[shift + 2*Nbpts_D*nbpts_linelets  + noind*nbpts_linelets : shift + 2*Nbpts_D*nbpts_linelets + (noind+1)*nbpts_linelets])

                                     shift = shift + Nbpts_D*(nbpts_linelets*nfields_lin + 1)


         FastC.HOOK['linelets_real'] = linelets

         Internal.createUniqueChild(tc, 'linelets_real', 'DataArray_t', FastC.HOOK['linelets_real'])
         Internal.createUniqueChild(tc, 'linelets_int', 'DataArray_t', FastC.HOOK['linelets_int'])

      return None

#==============================================================================
# Compute linelets info
# Adaptive h0
#==============================================================================
def computeLineletsInfo(tc, Re=6.e6, Cf_law='ANSYS', Lref=1., q=1.2):
    if Cf_law == 'ANSYS':
        def compute_Cf(Re):
            return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        def compute_Cf(Re):
            return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        def compute_Cf(Re):
            return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        def compute_Cf(Re):
            return 1.46*Re**(-0.5)

    hmod = 0.

    zones_tc = Internal.getZones(tc)

    for z in zones_tc:
        subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        if subRegions is not None:
              for s in subRegions:
                zsrname = s[0]
                sname  = zsrname[0:2]
                if sname == 'IB':
                  zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                  yI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                  xI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                  zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                  yW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                  xW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                  pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                  Nbpts_D       = numpy.shape(pointlistD[1])[0]

                  for i in range(Nbpts_D):
                    hmod = max(hmod, math.sqrt( (xI[i]-xW[i])**2 + (yI[i]-yW[i])**2 + (zI[i]-zW[i])**2 ))


    h0 = (Lref*math.sqrt(2))/(Re*math.sqrt(compute_Cf(Re)))
    # hmod = sum(0,n) h0*q**n
    n = int(math.ceil(math.log(1-((hmod/h0)*(1-q)))/(math.log(q)))) - 1
    # (n+1) intervalles / (n+2) points

    return (h0, h0*q**n, n+2)   

#==============================================================================
# Compute linelets info
# Fixed h0 1e-6
#==============================================================================
def computeLineletsInfo2(tc, q=1.2):
    hmod = 0.

    zones_tc = Internal.getZones(tc)

    for z in zones_tc:
        subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        if subRegions is not None:
              for s in subRegions:
                zsrname = s[0]
                sname  = zsrname[0:2]
                if sname == 'IB':
                  zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')[1]
                  yI = Internal.getNodeFromName1(s,'CoordinateY_PI')[1]
                  xI = Internal.getNodeFromName1(s,'CoordinateX_PI')[1]
                  zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')[1]
                  yW = Internal.getNodeFromName1(s,'CoordinateY_PW')[1]
                  xW = Internal.getNodeFromName1(s,'CoordinateX_PW')[1]
                  pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                  Nbpts_D       = numpy.shape(pointlistD[1])[0]

                  for i in range(Nbpts_D):
                    hmod = max(hmod, math.sqrt( (xI[i]-xW[i])**2 + (yI[i]-yW[i])**2 + (zI[i]-zW[i])**2 ))
                    
    h0 = 1.e-6
    n = int(math.ceil(math.log(1-((hmod/h0)*(1-q)))/(math.log(q)))) - 1

    return (h0, h0*q**n, n+2)   


#==============================================================================
# Prepare IBC ODE (Spalart 1D)
# Info stored and compacted in IBCD ZONES
# Para friendly
#==============================================================================
def _createTBLESA2(t, tc, h0, hn, nbpts_linelets=45):
      nfields_lin  = 6 # u,nutild,y,matm,mat,matp
      count_racIBC = 0
      addr_IBC = 0
      addrIBC  = []
      size_IBC = 0
      isODEDataPresent = 0
      distIW = 0.

      zones_tc = Internal.getZones(tc)

      for z in zones_tc:
          subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
          if subRegions is not None:
                for s in subRegions:
                  zsrname = s[0]
                  sname  = zsrname[0:2]
                  if sname == 'IB':
                     zsrname = zsrname.split('_')
                     if len(zsrname)<3:
                        print('Warning: createTBLESA: non consistent with the version of IBM preprocessing.')
                     else:
                       if zsrname[1] == '16':
                        print('Using MuskerSA2 wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1
                       elif zsrname[1] == '17':
                        print('Using FULL_TBLE_SA wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1
                       elif zsrname[1] == '18':
                        print('Using FULL_TBLE_Prandtl wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1
                       elif zsrname[1] == '19':
                        print('Using MafzalSA wall model')
                        pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                        Nbpts_D       = numpy.shape(pointlistD[1])[0]
                        addrIBC.append(size_IBC)
                        size_IBC      = size_IBC + Nbpts_D*(nbpts_linelets*nfields_lin + 1)
                        count_racIBC  = count_racIBC + 1
                        u_ODE = Internal.getNodeFromName1(s, 'VelocityT_ODE')
                        if u_ODE is not None: isODEDataPresent = 1

      if size_IBC > 0:
        if not isODEDataPresent:   # No ODE Data in tc
         print('ODE data missing in tc, they will be built...')

         if Internal.getNodeFromName1(t, 'nbpts_linelets') is None: 
          Internal.createUniqueChild(t, 'nbpts_linelets', 'DataArray_t', value=nbpts_linelets)

         # Nbpts_Dtot = size_IBC/(nbpts_linelets*nfields_lin)
         Nbpts_Dtot = size_IBC//(nbpts_linelets*nfields_lin+1)
         _addrIBC = numpy.asarray([nbpts_linelets,count_racIBC-1,Nbpts_Dtot] + addrIBC + [0]*Nbpts_Dtot ,dtype=numpy.int32)

         linelets = numpy.zeros(size_IBC + Nbpts_Dtot, dtype=numpy.float64)

         shift = 0
         for z in zones_tc:
           subRegions =  Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
           if subRegions is not None:
              for s in subRegions:
                  zsrname = s[0]
                  sname  = zsrname[0:2]
                  if sname == 'IB':
                     zsrname = zsrname.split('_')
                     if len(zsrname)<3:
                       print('Warning: createTBLESA2: non consistent with the version of IBM preprocessing.')
                     else:
                       if zsrname[1] == '16' or zsrname[1] == '17' or zsrname[1] == '18' or zsrname[1] == '19':
                         pointlistD    = Internal.getNodeFromName1(s, 'PointListDonor')
                         Nbpts_D       = numpy.shape(pointlistD[1])[0]

                         gradxP = Internal.getNodeFromName1(s, 'gradxPressure')
                         if gradxP is None:
                            Internal.createUniqueChild(s, 'gradxPressure', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))
                            Internal.createUniqueChild(s, 'gradyPressure', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))
                            Internal.createUniqueChild(s, 'gradzPressure', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))


                         zI = Internal.getNodeFromName1(s,'CoordinateZ_PI')
                         valzI = Internal.getValue(zI)
                         yI = Internal.getNodeFromName1(s,'CoordinateY_PI')
                         valyI = Internal.getValue(yI)
                         xI = Internal.getNodeFromName1(s,'CoordinateX_PI')
                         valxI = Internal.getValue(xI)
                         zW = Internal.getNodeFromName1(s,'CoordinateZ_PW')
                         valzW = Internal.getValue(zW)
                         yW = Internal.getNodeFromName1(s,'CoordinateY_PW')
                         valyW = Internal.getValue(yW)
                         xW = Internal.getNodeFromName1(s,'CoordinateX_PW')
                         valxW = Internal.getValue(xW)

                         if Nbpts_D > 1:

                             for noind in range(0,Nbpts_D):

                                 b0 = valxI[noind]-valxW[noind]
                                 b1 = valyI[noind]-valyW[noind]
                                 b2 = valzI[noind]-valzW[noind]
                                 normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)

                                 if hn < 0:
                                   hn = 0.6*normb

                                 # _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                 _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets], nbpts_linelets, 0.0, normb, h0, hn, 2)
                         else:

                                 b0 = valxI-valxW
                                 b1 = valyI-valyW
                                 b2 = valzI-valzW
                                 normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)

                                 if hn < 0:
                                   hn = 0.6*normb

                                 # _stretch(linelets[shift : shift + nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                 _stretch(linelets[shift : shift + nbpts_linelets], nbpts_linelets, 0.0, normb, h0, hn, 2)


                         Internal.createUniqueChild(s, 'CoordinateN_ODE', 'DataArray_t',
                                                    linelets[shift : shift + Nbpts_D*nbpts_linelets])
                         Internal.createUniqueChild(s, 'VelocityT_ODE', 'DataArray_t',
                                                    linelets[shift + 1*Nbpts_D*nbpts_linelets :  shift + 2*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'TurbulentSANuTilde_ODE', 'DataArray_t',
                                                    linelets[shift + 2*Nbpts_D*nbpts_linelets :  shift + 3*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'Psi_ODE', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D*nbpts_linelets, numpy.float64))
                         Internal.createUniqueChild(s, 'Matm_ODE', 'DataArray_t',
                                                    linelets[shift + 3*Nbpts_D*nbpts_linelets :  shift + 4*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'Mat_ODE', 'DataArray_t',
                                                    linelets[shift + 4*Nbpts_D*nbpts_linelets :  shift + 5*Nbpts_D*nbpts_linelets ])
                         Internal.createUniqueChild(s, 'Matp_ODE', 'DataArray_t',
                                                    linelets[shift + 5*Nbpts_D*nbpts_linelets :  shift + 6*Nbpts_D*nbpts_linelets ])

                         Internal.createUniqueChild(s, 'alphasbeta_ODE', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))
                         Internal.createUniqueChild(s, 'index_ODE', 'DataArray_t',
                                                    numpy.zeros(Nbpts_D, numpy.float64))

                         # Internal.createUniqueChild(s, 'Matm_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D*45, numpy.float64)*2)

                         # Internal.createUniqueChild(s, 'Mat_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D*45, numpy.float64)*3)
                         
                         # Internal.createUniqueChild(s, 'Matp_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D*45, numpy.float64)*4)

                         # Internal.createUniqueChild(s, 'alphasbeta_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D, numpy.float64)*5)

                         # Internal.createUniqueChild(s, 'index_ODE', 'DataArray_t',
                         #                            numpy.ones(Nbpts_D, numpy.float64)*6)


                         shift = shift + Nbpts_D*(nbpts_linelets*nfields_lin + 1)

      return None



#==============================================================================
# Update gradx/y/zPressure in tc IBCD zones using first (slow) method without compact
# For tc2 the type may have to be '2_IBCD' instead of 'IBCD'
#==============================================================================
def _updateGradPInfo(t, tc, metrics, type='IBCD'): 

  varGrad = ['Density', 'Temperature']
  variablesIBC = ["Density", "Temperature", "gradxDensity", "gradyDensity", "gradzDensity", "gradxTemperature", "gradyTemperature", "gradzTemperature"]

  for var in varGrad:
   C._rmVars(t, 'gradx{}'.format(var))
   C._rmVars(t, 'grady{}'.format(var))
   C._rmVars(t, 'gradz{}'.format(var))

   _computeGrad(t, metrics, varGrad)
   for v in variablesIBC: C._cpVars(t, 'centers:'+v, tc, v)
   XOD._setInterpTransfers(t, tc, variables=variablesIBC, variablesIBC=[], compact=0, compactD=0)
   XOD._setIBCTransfers4GradP(t, tc, variablesIBC=variablesIBC)

   C._rmVars(tc, variablesIBC)

  return None

#==============================================================================
# Update gradx/y/zPressure in tc IBCD zones using first (slow) method without compact
# For tc2 the type may have to be '2_IBCD' instead of 'IBCD'
#==============================================================================
def _updateGradPInfoHO(t, tc, metrics, type='IBCD', algo="Fast"): 
  import Post.PyTree as P

  variablesIBC = ["Pressure"]
  varGrad = []

  for var in variablesIBC:
    for direction in ["x", "y", "z"]:
      C._rmVars(t, "grad{}{}".format(direction,var))
      varGrad.append("grad{}{}".format(direction,var))

  varGrad2 = []
  for var in varGrad:
    for direction in ["x", "y", "z"]:
      C._rmVars(t, "grad{}{}".format(direction,var))
      varGrad2.append("grad{}{}".format(direction,var))

  if algo == "Fast":
    _computeGrad(t, metrics, variablesIBC)
    _computeGrad(t, metrics, varGrad)
  else:
    for var in variablesIBC:
      print(var)
      t = P.computeGrad(t, 'centers:'+var)
    for var in varGrad:
      print(var)
      t = P.computeGrad(t, 'centers:'+var)

  variables = variablesIBC+varGrad+varGrad2

  print(variables)

  for v in variables: C._cpVars(t, 'centers:'+v, tc, v)
  XOD._setInterpTransfers(t, tc, variables=variables, variablesIBC=[], compact=0, compactD=0)
  XOD._setIBCTransfers4GradP2(t, tc, variablesIBC=variables)

  C._rmVars(tc, variables)

  return None

def _UpdateUnsteadyJoinParam(t, tc, omega, timelevelInfos, split='single', root_steady='tc_steady', root_unsteady='tc_', dir_steady='.', dir_unsteady='.', init=False):

    #on cree les noeud infos insta pour chimere perio s'il n'existe pas 
    TimeLevelOpts=['TimeLevelMotion','TimeLevelTarget','Iteration'] 
    for opt in TimeLevelOpts:
       tmp = Internal.getNodeFromName1(t, opt)
       if tmp is None: Internal.createUniqueChild(t, opt, 'DataArray_t', value=0)

    #mise a jour des info temporelle
    if not init:
      Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0] +=1
      Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0] +=1
      Internal.getNodeFromName1(t, 'Iteration'      )[1][0] +=1
 
    timelevel_motion = Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0]
    timelevel_target = Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0]
    iteration        = Internal.getNodeFromName1(t, 'Iteration')[1][0]

    timelevel_period = timelevelInfos["TimeLevelPeriod"]
    timelevel_360    = timelevelInfos["TimeLevel360"]
    timelevel_perfile= timelevelInfos["TimeLevelPerFile"]
    timelevel_axeRot = timelevelInfos["TimeLevelRotationAxis"]

    No_period = timelevel_motion//timelevel_period 
    #
    #target in no more in tc; need need data in a new file
    #
    #if timelevel_target == timelevel_perfile+1 or tc is None:
    if timelevel_target == timelevel_perfile or tc is None:

       tmp  = No_period*timelevel_period
       root = timelevel_perfile + ( (timelevel_motion - tmp)//timelevel_perfile)*timelevel_perfile

       # steady part
       FILE = dir_steady +'/'+ root_steady + '.cgns'
       if os.access(FILE, os.F_OK): tc  = C.convertFile2PyTree(FILE)
       else: print("error reading %s."%FILE)

       # unsteady part
       FILE = dir_unsteady +'/'+root_unsteady+str(root)+'.cgns'
       if os.access(FILE, os.F_OK):
           tc_inst  = C.convertFile2PyTree(FILE)
           print('File inst=', FILE, 'target=', timelevel_target, 'motionlevel=', timelevel_motion)
       else: print("error reading %s."%FILE)

       #mise a zero timelevel_target en fin de fichier ou fin de periode azymutale
       if timelevel_target == timelevel_perfile or timelevel_motion%timelevel_period == 0:
          timelevel_target =0
          Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0] = timelevel_target

       #
       #timelevel_motion larger than calculated peridicity; need to modify angle of rotation for azymuth periodicity
       #
       if timelevel_motion >= timelevel_period:

           No_period     = timelevel_motion//timelevel_period
           iteration_loc = timelevel_motion - No_period*timelevel_period 

           bases  = Internal.getNodesFromType1(tc_inst , 'CGNSBase_t')       # noeud
           Rotors  = ['Base02', 'Base04', 'Base06']
           Stators = ['Base01', 'Base03', 'Base05', 'Base07']

           sign =-1
           if omega > 0: sign = 1
           for base in bases:
             if   base[0] in Rotors  : teta = -2*math.pi*timelevel_period/timelevel_360*No_period*sign
             elif base[0] in Stators : teta =  2*math.pi*timelevel_period/timelevel_360*No_period*sign
             zones  = Internal.getNodesFromType1(base , 'Zone_t')       # noeud
             for z in zones:
               angles = Internal.getNodesFromName2(z, 'RotationAngle')
               for angle in angles: angle[1][:]= angle[1][:] + teta*timelevel_axeRot[:]
       else:
            iteration_loc = timelevel_motion

       tc = Internal.merge( [tc, tc_inst] )

       # Get omp_mode
       ompmode = FastC.OMP_MODE
       node = Internal.getNodeFromName1(t, '.Solver#define')
       if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

       # Reordone les zones pour garantir meme ordre entre t et tc
       FastC._reorder(t, tc, ompmode)

       # Compactage arbre transfert
       g = None; l = None
       zones = Internal.getZones(t)

       X.miseAPlatDonorTree__(zones, tc, graph=g, list_graph=l)

    #
    #timelevel_motion larger than number of timelevels for 360degre
    #
    if timelevel_motion > timelevel_360:
       timelevel_motion = 0
       Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0] =timelevel_motion

    own  = Internal.getNodeFromName1(t, '.Solver#ownData')
    if own is not None: 
       dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')[1]
       dtloc[3] = timelevel_motion
       dtloc[4] = timelevel_target

    return tc, graph

#==============================================================================
def checkBalance(t):
    zones = Internal.getZones(t)

    size_reelle=0
    size_transf=0
    for z in zones:
       dim = Internal.getZoneDim(z)
       iv  = dim[1]-5
       jv  = dim[2]-5
       kv  = dim[3]-5
       trans =  4*(iv*jv + iv*kv + kv*jv)
       size_reelle = size_reelle +iv*jv*kv
       size_transf = size_transf + trans
    print("size Pb: totale =", size_reelle+size_transf, 'woghost=', size_reelle, 'NB ghost=',size_transf, 'Nbzone=',len(zones))
    return None

#==============================================================================
# Interface for Vtune/Advisor collection control
#==============================================================================
def itt(var):
    if var == 'pause': ivar =1
    else: ivar = 0
    print("itt collection (Vtune/Advisor)", var)
    fasts.itt(ivar)
    return None

#==============================================================================
def _applyBC(t, metrics, hook1, nstep, var="Density"):
    zones = Internal.getZones(t)

    fasts._applyBC(zones, metrics, hook1, nstep, var)
    return None

#==============================================================================
def _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, omp_mode, hook1, nitmax=1, rk=1, exploc=0, num_passage=1, gradP=False, TBLE=False, isWireModel=False):
   # timecount = numpy.zeros(4, dtype=numpy.float64)
   
   varsGrad = []
   if hook1['lexit_lu'] ==0:
 
       #transfert
       if tc is not None:
           tc_compact = Internal.getNodeFromName1( tc, 'Parameter_real')
           #Si param_real n'existe pas, alors pas de raccord dans tc
           if tc_compact is not  None:

              param_real= tc_compact[1]
              param_int = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1["neq_max"] == 5: varType = 2
              else                    : varType = 21

              dtloc = hook1['dtloc']

              for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              if gradP:
                  varsGrad = ['Density', 'gradxDensity']
                  for v in varsGrad: C._cpVars(zones, 'centers:'+v, zonesD,  v)
                  varType = 22 ; type_transfert = 2 ; no_transfert   = 1
                  # varType 22 means that there is 6 variables and 6 gradVariables to transfers
                  # to avoid interfering with the original ___setInterpTransfers, we instead use ___setInterpTransfers4GradP
                  Connector.connector.___setInterpTransfers4GradP(zones, zonesD, varsGrad, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage)#,timecount)
                  varType = 21

              if TBLE: # need to be optimised
                  variablesIBC = ["Density", "Temperature", "gradxDensity", "gradyDensity", "gradzDensity", "gradxTemperature", "gradyTemperature", "gradzTemperature",
                  'gradxVelocityX','gradyVelocityX','gradzVelocityX',
                  'gradxVelocityY','gradyVelocityY','gradzVelocityY',
                  'gradxVelocityZ','gradyVelocityZ','gradzVelocityZ',
                  'VelocityX','VelocityY','VelocityZ',
                  ]
                  for v in variablesIBC: C._cpVars(zones, 'centers:'+v, zonesD, v)
                  XOD._setInterpTransfers(zones, zonesD, type_transfert=0, variables=variablesIBC, variablesIBC=[], compact=0)
                  XOD._setIBCTransfers4FULLTBLE(zones, zonesD, variablesIBC=variablesIBC)

              
              if isWireModel:
                  type_transfert  = 1
                  no_transfert    = 1
                  isWireModel_int = 1
                  Connector.connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage, isWireModel_int)#,timecount)
                  isWireModel_int = 2
                  Connector.connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage, isWireModel_int)#,timecount)

              type_transfert  = 2  # 0= ID uniquement, 1= IBC uniquement, 2= All
              no_transfert    = 1  # dans la list des transfert point a point
              isWireModel_int = 0
              Connector.connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage, isWireModel_int)#,timecount)
                  
       #apply BC
       #t0=timeit.default_timer()
       if exploc != 1:
          _applyBC(zones, metrics, hook1, nstep, var=vars[0])
       #t1=timeit.default_timer()
       #print("Time BC",(t1-t0))

   return None

#==============================================================================
# modified _fillGhostCells for two image points (tc + tc2) for gradP 
#==============================================================================
def _fillGhostcells2(zones, tc, tc2, metrics, timelevel_target, vars, nstep, omp_mode, hook1, nitmax=1, rk=1, exploc=0, num_passage=1, gradP=False, TBLE=False):

   # timecount = numpy.zeros(4, dtype=numpy.float64)
   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None:
           tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
           #Si param_real n'existe pas, alors pas de raccord dans tc
           if tc_compact is not None:

              param_real= tc_compact[1]
              param_int = Internal.getNodeFromName1(tc, 'Parameter_int')[1]
              zonesD    = Internal.getZones(tc)
              zonesD2   = Internal.getZones(tc2)

              if hook1["neq_max"] == 5: varType = 2
              else                    : varType = 21

              dtloc = hook1['dtloc']

              for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              if gradP:
                  tc2_compact = Internal.getNodeFromName1( tc2, 'Parameter_real')
                  param_real2= tc2_compact[1]
                  param_int2 = Internal.getNodeFromName1(tc2, 'Parameter_int' )[1]
                  varsGrad = ['Density', 'gradxDensity']
                  for v in varsGrad: 
                    C._cpVars(zones, 'centers:'+v, zonesD,  v)
                    C._cpVars(zones, 'centers:'+v, zonesD2, v)
                  #tc2 -> RCV ZONES
                  varType = 23 ; type_transfert = 2 ; no_transfert   = 1
                  Connector.connector.___setInterpTransfers4GradP(zones, zonesD2, varsGrad, param_int2, param_real2, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage)
                  #RCV ZONES -> tc
                  varType = 24 ; type_transfert = 1 ; no_transfert   = 1
                  Connector.connector.___setInterpTransfers4GradP(zones, zonesD, varsGrad, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage)
                  varType = 21

              elif TBLE: # need to be optimised
                  variablesIBC = ["Density", "Temperature", "gradxDensity", "gradyDensity", "gradzDensity", "gradxTemperature", "gradyTemperature", "gradzTemperature",
                  'gradxVelocityX','gradyVelocityX','gradzVelocityX',
                  'gradxVelocityY','gradyVelocityY','gradzVelocityY',
                  'gradxVelocityZ','gradyVelocityZ','gradzVelocityZ',
                  'VelocityX','VelocityY','VelocityZ']
                  for v in variablesIBC: C._cpVars(zones, 'centers:'+v, zonesD, v)
                  XOD._setInterpTransfers(zones, zonesD,  type_transfert=0, variables=variablesIBC, variablesIBC=[], compact=0)
                  for v in variablesIBC: C._cpVars(zones, 'centers:'+v, zonesD2, v)
                  XOD._setIBCTransfers4FULLTBLE(zones, zonesD2, variablesIBC=variablesIBC)                  
                  XOD._setIBCTransfers4FULLTBLE2(zones, zonesD, variablesIBC=variablesIBC)                  


              type_transfert = 2  # 0= ID uniquement, 1= IBC uniquement, 2= All
              no_transfert   = 1  # dans la list des transfert point a point
              Connector.connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage, 0)#,timecount)


       #apply BC
       #t0=timeit.default_timer()
       if exploc != 1:
       #if rk != 3 and exploc != 2:
          _applyBC(zones, metrics, hook1, nstep, var=vars[0])
       #t1=timeit.default_timer()
       #print("Time BC",(t1-t0))

   return None

#==============================================================================
# Cree un noeud POST
# IN: t: tree
# IN: dir: direction de la moyenne '0', 'i', 'j', 'k', 'ij', 'ik', 'jk'
# IN: vars: variables concernees
# IN: nsample: nbre d'echantillons dans la moyenne
# OUT: return a tree with POST node
#==============================================================================
def createStatNodes(t, dir='0', vars=[], nsamples=0):
    """Create node in tree to store stats."""
    try: import Transform.PyTree as T
    except: raise ImportError("createStatNodes: requires Transform module.")

    PostBaseName = 'POST' # nom de la base POST
    DataNodeName = '.Solver#post'
    vars0 = ['CoordinateX','CoordinateY','CoordinateZ']

    tmy = C.newPyTree([PostBaseName])

    b = Internal.getNodesFromName1(tmy, PostBaseName)

    varmy = ['MomentumX','MomentumY','MomentumZ','Density','Pressure','Pressure^2','ViscosityEddy','rou^2','rov^2','row^2','rouv','rouw','rovw']
    lgrad = 0
    for var in vars:
      if var == 'cylindrique' or var == 'cylx':
        varmy[1] = 'Momentum_t'
        varmy[2] = 'Momentum_r'
        varmy[8] = 'roU_t^2'
        varmy[9] = 'roU_r^2'
        varmy[10] = 'rouU_t'
        varmy[11] = 'rouU_r'
        varmy[12] = 'roU_tU_r'
      elif var == 'cylz':
        varmy[0] = 'Momentum_t'
        varmy[1] = 'Momentum_r'
        varmy[7] = 'roU_t^2'
        varmy[8] = 'roU_r^2'
        varmy[10] = 'roU_tU_r'
        varmy[11] = 'rowU_t'
        varmy[12] = 'rowU_r'
      elif var == 'thermique':
        varmy += ['Temperature','T^2','rouT','rovT','rowT','Eps_T']
        lgrad =  1

    for i in range(len(varmy)): varmy[i] = 'centers:'+varmy[i]

    # on determine le nbr de cellule fictive active pour le calcul des moyennes
    numcellfic = 2
    ific       = 2   # a adapter en DF
    if lgrad == 1: numcellfic = 1

    zones = []
    for b0 in Internal.getNodesFromType1(t,'CGNSBase_t'):
        if b0[0] != PostBaseName:
            zones += Internal.getZones(b0)

    own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    tmy[2].append(own)

    if dir == '0':
        for z in zones:
            #
            datap = numpy.empty((17), numpy.int32)
            #
            dim  = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0
            zp = T.subzone(z, (1,1,1), (-1,-1,-1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)

            dim_my = Internal.getZoneDim(zp)

            datap[0]  = 0                                                # nbr direction homogene
            datap[1]  = 0                                                # dir homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1]-1)*(dim_my[2]-1)*(dim_my[3]-1)        # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = 1                                                # nbre cellule homogene
            datap[6]  = 1                                                # incerment i
            datap[7]  = dim_my[1]-1                                      # increment j
            datap[8]  = (dim_my[1]-1)*(dim_my[2]-1)                       # increment k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim_my[3] == 2: datap[10] =  0
            datap[11] =  1 - numcellfic                                  # loop moyenne imin
            datap[12] =  (dim_my[1]-1) -2*ific + numcellfic              # loop moyenne imax
            datap[13] =  1 - numcellfic                                  # loop moyenne jmin
            datap[14] =  dim_my[2]-1 -2*ific + numcellfic                   # loop moyenne jmax
            datap[15] =  1 - numcellfic*inck                             # loop moyenne kmin
            datap[16] =  dim_my[3]-1 -2*ific*inck + numcellfic*inck         # loop moyenne kmax
            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
            param_int[1][66]= 0                           # pas de padding pour les variables stats

            b[0][2].append(zp)

    elif dir == 'i':
        for z in zones:
            datap = numpy.empty((17), numpy.int32)
            dim   = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0
            zp = T.subzone(z, (1,1,1), (1,-1,-1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)

            dim_tr = Internal.getZoneDim(zp)
            dim_my = [ 0, 1 , dim_tr[1], dim_tr[2] ]
            datap[0]  = 1                                                # nbr direction homogene
            datap[1]  = 1                                                # dir homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1])*(dim_my[2]-1)*(dim_my[3]-1)        # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = dim[1]-1 -2*ific                                 # nbre cellule homogene
            datap[6]  = 0                                                # incerment i
            datap[7]  = dim_my[1]-1                                      # increment j
            datap[8]  =(dim_my[1]-1)*(dim_my[2]-1)                       # increment k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim[3] == 2: datap[10] =  0
            datap[11] =  1                                               # loop moyenne imin
            datap[12] =  1                                               # loop moyenne imax
            datap[13] =  1 - numcellfic                                  # loop moyenne jmin
            datap[14] =  datap[7] -2*ific + numcellfic                   # loop moyenne jmax
            datap[15] =  1 - numcellfic*inck                             # loop moyenne kmin
            datap[16] =  datap[8] -2*ific*inck + numcellfic*inck         # loop moyenne kmax

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS

            param_int[1][0] = dim_my[1]             #nijk
            param_int[1][1] = dim_my[2]-1
            param_int[1][2] = dim_my[3]-1
            param_int[1][3] = datap[9]
            param_int[1][4] = datap[10]
            param_int[1][20]= 1                     #ijkv
            param_int[1][21]= dim_my[2]-1 -2*ific
            param_int[1][22]= dim_my[3]-1 -2*ific*inck
            param_int[1][66]= 0                           # pas de padding pour les variables stats

            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])
            b[0][2].append(zp)

    elif dir == 'j':
        for z in zones:
            #
            datap = numpy.empty((17), numpy.int32)
            #
            dim  = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0

            zp = T.subzone(z, (1,1,1), (-1,1,-1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)

            dim_tr = Internal.getZoneDim(zp)
            dim_my = [ 0, dim_tr[1], 1, dim_tr[2] ]
            datap[0]  = 1                                                # nbr direction homogene
            datap[1]  = 2                                                # dir homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1]-1)*(dim_my[2])*(dim_my[3]-1)        # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = dim[2]-1  -2*ific                                # nbre cellule homogene
            datap[6]  = 1                                                # incerment i
            datap[7]  = 0                                                # increment j
            datap[8]  =(dim_my[1]-1)                                     # increment k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim[3] == 2: datap[10] =  0
            datap[11] =  1 - numcellfic                                  # loop moyenne imin
            datap[12] =  datap[6] -2*ific + numcellfic                   # loop moyenne imax
            datap[13] =  1                                               # loop moyenne jmin
            datap[14] =  1                                               # loop moyenne jmax
            datap[15] =  1 - numcellfic*inck                             # loop moyenne kmin
            datap[16] =  datap[8] -2*ific*inck + numcellfic*inck         # loop moyenne kmax

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
            param_int[1][0] = dim_my[1]-1             #nijk
            param_int[1][1] = dim_my[2]
            param_int[1][2] = dim_my[3]-1
            param_int[1][3] = datap[9]
            param_int[1][4] = datap[10]
            param_int[1][20]= dim_my[2]-1 -2*ific   #ijkv
            param_int[1][21]= 1
            param_int[1][22]= dim_my[3]-1 -2*ific*inck
            param_int[1][66]= 0                           # pas de padding pour les variables stats

            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])
            b[0][2].append(zp)

    elif dir == 'k':
        for z in zones:
            datap = numpy.empty(17, numpy.int32)

            dim  = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0

            zp = T.subzone(z, (1,1,1), (-1,-1,1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)
            dim_tr = Internal.getZoneDim(zp)
            dim_my = [ 0, dim_tr[1], dim_tr[2], 1 ]

            datap[0]  = 1                                                # nbr direction homogene
            datap[1]  = 3                                                # dir homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1]-1)*(dim_my[2]-1)*(dim_my[3])          # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = dim[3]-1  -2*ific*inck                           # nbre cellule homogene
            datap[6]  = 1                                                # nbre cellule i pour calcul adress
            datap[7]  = dim_my[1]-1                                      # nbre cellule j
            datap[8]  = 0                                                # nbre cellule k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim[3] == 2: datap[10] =  0
            datap[11] =  1 - numcellfic                                  # loop moyenne imin
            datap[12] =  (dim_my[1]-1) -2*ific + numcellfic              # loop moyenne imax
            datap[13] =  1 - numcellfic                                  # loop moyenne jmin
            datap[14] =  (dim_my[2]-1) -2*ific + numcellfic              # loop moyenne jmax
            datap[15] =  1                                               # loop moyenne kmin
            datap[16] =  1                                               # loop moyenne kmax

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
            param_int[1][0] = dim_my[1]-1           #nijk
            param_int[1][1] = dim_my[2]-1
            param_int[1][2] = dim_my[3]
            param_int[1][3] = datap[9]
            param_int[1][4] = 0
            param_int[1][20]= dim_my[1]-1 -2*ific   #ijkv
            param_int[1][21]= dim_my[2]-1 -2*ific
            param_int[1][22]= 1
            param_int[1][66]= 0                           # pas de padding pour les variables stats

            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])
            b[0][2].append(zp)

    elif dir == 'ij':
        for z in zones:

            datap = numpy.empty((17), numpy.int32)
            dim   = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0

            zp = T.subzone(z, (1,1,1), (1,1,-1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)

            dim_tr = Internal.getZoneDim(zp)
            dim_my = [ 0, 1, 1, dim_tr[1] ]
            datap[0]  = 2                                                # nbr direction homogene
            datap[1]  = 3                                                # dir non homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1])*(dim_my[2])*(dim_my[3]-1)        # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = (dim[1]-1  -2*ific)*(dim[2]-1  -2*ific)          # nbre cellule homogene
            datap[6]  = 0                                                # nbre cellule i pour calcul adress
            datap[7]  = 0                                                # nbre cellule j
            datap[8]  = 1                                                # nbre cellule k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim[3] == 2: datap[10] =  0
            datap[11] =  1                                               # loop moyenne imin
            datap[12] =  1                                               # loop moyenne imax
            datap[13] =  1                                               # loop moyenne jmin
            datap[14] =  1                                               # loop moyenne jmax
            datap[15] =  1 - numcellfic*inck                             # loop moyenne kmin
            datap[16] =  datap[8] -2*ific*inck + numcellfic*inck         # loop moyenne kmax

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
            param_int[1][0] = dim_my[1]               #nijk
            param_int[1][1] = dim_my[2]
            param_int[1][2] = dim_my[3]-1
            param_int[1][3] = datap[9]
            param_int[1][4] = datap[10]
            param_int[1][20]= 1                     #ijkv
            param_int[1][21]= 1
            param_int[1][22]= dim_my[3]-1 -2*ific*inck
            param_int[1][66]= 0                           # pas de padding pour les variables stats

            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])
            b[0][2].append(zp)

    elif dir == 'ik':
        for z in zones:
            datap = numpy.empty((17), numpy.int32)

            dim   = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0

            zp = T.subzone(z, (1,1,1), (1,-1,1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)

            dim_tr = Internal.getZoneDim(zp)
            dim_my = [ 0, 1, dim_tr[1], 1 ]
            datap[0]  = 2                                                # nbr direction homogene
            datap[1]  = 2                                                # dir non homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1])*(dim_my[2]-1)*(dim_my[3])        # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = (dim[1]-1  -2*ific)*(dim[3]-1  -2*ific*inck)     # nbre cellule homogene
            datap[6]  = 0                                                # nbre cellule i pour calcul adress
            datap[7]  = 1                                                # nbre cellule j
            datap[8]  = 0                                                # nbre cellule k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim[3] == 2: datap[10] =  0
            datap[11] =  1                                               # loop moyenne imin
            datap[12] =  1                                               # loop moyenne imax
            datap[13] =  1 - numcellfic                                  # loop moyenne jmin
            datap[14] =  datap[7] -2*ific + numcellfic                   # loop moyenne jmax
            datap[15] =  1                                               # loop moyenne kmin
            datap[16] =  1                                               # loop moyenne kmax

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
            param_int[1][0] = dim_my[1]               #nijk
            param_int[1][1] = dim_my[2]-1
            param_int[1][2] = dim_my[3]
            param_int[1][3] = datap[9]
            param_int[1][4] = datap[10]
            param_int[1][20]= 1                     #ijkv
            param_int[1][21]= dim_my[2]-1 -2*ific
            param_int[1][22]= 1
            param_int[1][66]= 0                           # pas de padding pour les variables stats

            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])
            b[0][2].append(zp)

    elif dir == 'jk':
        for z in zones:
            datap = numpy.empty((17), numpy.int32)
            dim   = Internal.getZoneDim(z)
            inck = 1
            if dim[3] == 2: inck =  0

            zp = T.subzone(z, (1,1,1), (-1,1,1)) ; zp[0] = z[0]
            C._extractVars(zp, vars0)
            for var in varmy: C._initVars(zp, var, 0.)

            dim_tr = Internal.getZoneDim(zp)
            dim_my = [ 0, dim_tr[1], 1, 1 ]
            datap[0]  = 2                                                # nbr direction homogene
            datap[1]  = 1                                                # dir non homogne
            datap[2]  = nsamples                                         # nbre echantillon
            datap[3]  = (dim_my[1]-1)*(dim_my[2])*(dim_my[3])            # nbre d'element par champ
            datap[4]  = len( varmy )                                     # nbre champ statistique
            datap[5]  = (dim[2]-1  -2*ific)*(dim[3]-1  -2*ific*inck)     # nbre cellule homogene
            datap[6]  = 1                                                # nbre cellule i pour calcul adress
            datap[7]  = 0                                                # nbre cellule j
            datap[8]  = 0                                                # nbre cellule k
            datap[9]  =  ific                                            # ific: a adapter  si DF
            datap[10] =  ific                                            # kfic: a adapter  si DF
            if dim_my[3] == 2: datap[10] =  0
            datap[11] =  1 - numcellfic                                  # loop moyenne jmin
            datap[12] =  datap[6] -2*ific + numcellfic                   # loop moyenne jmax
            datap[13] =  1                                               # loop moyenne imin
            datap[14] =  1                                               # loop moyenne imax
            datap[15] =  1                                               # loop moyenne kmin
            datap[16] =  1                                               # loop moyenne kmax

            param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%z[0])

            ##modifie le moeud param_int de l'arbre tmy (issue de subzone) pour la fonction initNuma
            param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
            param_int[1][0] = dim_my[1]-1             #nijk
            param_int[1][1] = dim_my[2]
            param_int[1][2] = dim_my[3]
            param_int[1][3] = datap[9]
            param_int[1][4] = datap[10]
            param_int[1][20]= dim_my[1]-1 -2*ific   #ijkv
            param_int[1][21]= 1
            param_int[1][22]= 1
            param_int[1][66]= 0                           # pas de padding pour les variables stats
            zp[2].append([DataNodeName,datap,[],'UserDefinedData_t'])
            b[0][2].append(zp)

    else: raise ValueError("createStatNodes: not a valid direction.")

    Internal._rmNodesByType(b,'ZoneBC_t')
    Internal._rmNodesByType(b,'ZoneGridConnectivity_t')

    FastC._compact(tmy, fields=varmy)

    return tmy

#==============================================================================
# compute statistique in place
#==============================================================================
def _computeStats(t, tmy, metrics):
    """Compute the space/time average of flowfields in tmy."""
    zones    = Internal.getZones(t)
    zones_my = Internal.getZones(tmy)

    node = Internal.getNodeFromName1(t, '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = FastC.OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    fasts.computePT_my(zones, zones_my, metrics, ompmode)
    return None

#==============================================================================
# compute statistique in place
#==============================================================================
def initStats(filename):

    tmy = C.convertFile2PyTree(filename)
    sol = Internal.getNodesFromName3(tmy, 'FlowSolution#Centers')
    var = Internal.getNodesFromType1(sol[0] , 'DataArray_t')

    varmy=[]
    for v in var: varmy.append('centers:'+v[0])
    FastC._compact(tmy, fields=varmy)

    return tmy

#==============================================================================
# compute statistique (U  et urms)
#==============================================================================
def _postStats(tmy):

    zones   = Internal.getZones(tmy)

    for z in zones:
      flowsol = Internal.getNodeFromName1(z,'FlowSolution#Centers')
      cyl = Internal.getNodeFromName1(flowsol,"Momentum_t")
      if cyl is None:
          vars=[['MomentumX','VelocityX'],['MomentumY','VelocityY'],['MomentumZ','VelocityZ']]
      else:
          vars=[['MomentumX','VelocityX'],['Momentum_t','Velocity_t'],['Momentum_r','Velocity_r']]

      #calcul vitesse
      dens= Internal.getNodeFromName1(flowsol,'Density')
      for v in vars:
         var = Internal.getNodeFromName1(flowsol,v[0])
         var[0]= v[1]
         var[1]= var[1]/dens[1]

      #calcul prms
      pres = Internal.getNodeFromName1(flowsol,'Pressure')
      pres2= Internal.getNodeFromName1(flowsol,'Pressure^2')
      pres2[0]='Prms'
      pres2[1]=numpy.sqrt( numpy.abs( pres2[1]-pres[1]**2))
      #calcul urms
      pres = Internal.getNodeFromName1(flowsol,'VelocityX')
      pres2= Internal.getNodeFromName1(flowsol,'rou^2')
      pres2[0]='Urms'
      pres2[1]=numpy.sqrt( numpy.abs(pres2[1]/dens[1] - pres[1]**2) )

      if cyl is None:
        #calcul vrms
        pres = Internal.getNodeFromName1(flowsol,'VelocityY')
        pres2= Internal.getNodeFromName1(flowsol,'rov^2')
        pres2[0]='Vrms'
        pres2[1]=numpy.sqrt( numpy.abs(pres2[1]/dens[1] - pres[1]**2) )
        #calcul wrms
        pres = Internal.getNodeFromName1(flowsol,'VelocityZ')
        pres2= Internal.getNodeFromName1(flowsol,'row^2')
        pres2[0]='Wrms'
        pres2[1]=numpy.sqrt( numpy.abs(pres2[1]/dens[1] - pres[1]**2) )
        #calcul rou'v'
        v1   = Internal.getNodeFromName1(flowsol,'VelocityX')
        v2   = Internal.getNodeFromName1(flowsol,'VelocityY')
        pres2= Internal.getNodeFromName1(flowsol,'rouv')
        pres2[0]="Rhou'v'"
        pres2[1]=pres2[1] - v1[1]*v2[1]*dens[1]
        #calcul rou'w'
        v1   = Internal.getNodeFromName1(flowsol,'VelocityX')
        v2   = Internal.getNodeFromName1(flowsol,'VelocityZ')
        pres2= Internal.getNodeFromName1(flowsol,'rouw')
        pres2[0]="Rhou'w'"
        pres2[1]=pres2[1] - v1[1]*v2[1]*dens[1]
        #calcul rov'w'
        v1   = Internal.getNodeFromName1(flowsol,'VelocityY')
        v2   = Internal.getNodeFromName1(flowsol,'VelocityZ')
        pres2= Internal.getNodeFromName1(flowsol,'rovw')
        pres2[0]="Rhov'w'"
        pres2[1]=pres2[1] - v1[1]*v2[1]*dens[1]

        C._initVars(z, '{centers:MomentumX}= {centers:Density}*{centers:VelocityX}')
        C._initVars(z, '{centers:MomentumY}= {centers:Density}*{centers:VelocityY}')
        C._initVars(z, '{centers:MomentumZ}= {centers:Density}*{centers:VelocityZ}')
      else:
        #calcul vrms
        pres = Internal.getNodeFromName1(flowsol,'Velocity_t')
        pres2= Internal.getNodeFromName1(flowsol,'roU_t^2')
        pres2[0]='Vt_rms'
        pres2[1]=numpy.sqrt( numpy.abs(pres2[1]/dens[1] - pres[1]**2) )
        #calcul wrms
        pres = Internal.getNodeFromName1(flowsol,'Velocity_r')
        pres2= Internal.getNodeFromName1(flowsol,'roU_r^2')
        pres2[0]='Vr_rms'
        pres2[1]=numpy.sqrt( numpy.abs(pres2[1]/dens[1] - pres[1]**2) )
        #calcul rou'v'
        v1   = Internal.getNodeFromName1(flowsol,'VelocityX')
        v2   = Internal.getNodeFromName1(flowsol,'Velocity_t')
        pres2= Internal.getNodeFromName1(flowsol,'rouU_t')
        pres2[0]="Rhou'v'"
        pres2[1]=pres2[1] - v1[1]*v2[1]*dens[1]
        #calcul rou'w'
        v1   = Internal.getNodeFromName1(flowsol,'VelocityX')
        v2   = Internal.getNodeFromName1(flowsol,'Velocity_r')
        pres2= Internal.getNodeFromName1(flowsol,'rouU_r')
        pres2[0]="Rhou'w'"
        pres2[1]=pres2[1] - v1[1]*v2[1]*dens[1]
        #calcul rov'w'
        v1   = Internal.getNodeFromName1(flowsol,'Velocity_t')
        v2   = Internal.getNodeFromName1(flowsol,'Velocity_r')
        pres2= Internal.getNodeFromName1(flowsol,'roU_tU_r')
        pres2[0]="Rhov'w'"
        pres2[1]=pres2[1] - v1[1]*v2[1]*dens[1]
        #calcul vrms
        C._initVars(z, '{centers:MomentumX}= {centers:Density}*{centers:VelocityX}')
        C._initVars(z, '{centers:MomentumY}= {centers:Density}*{centers:Velocity_t}')
        C._initVars(z, '{centers:MomentumZ}= {centers:Density}*{centers:Velocity_r}')

    return None

# Compact stat tree
def _compactStats(ts):
    sol = Internal.getNodesFromName3(ts, 'FlowSolution#Centers')
    var = Internal.getNodesFromType1(sol[0], 'DataArray_t')
    varmy = []
    for v in var: varmy.append('centers:'+v[0])
    FastC._compact(ts, fields=varmy)
    return None

# Return phase number
def phase(time, omega, dpsi):
    psi = time*omega*180./numpy.pi
    psi = psi %360
    p = int(psi/dpsi)
    return p

# IN: time: current time
# IN: omega: rotation speed (rad/s)
# IN: dpsi: phase range in psi
def _computePhaseStats(t, ts, metrics, time, omega, dpsi):
    import Converter.Mpi as Cmpi
    import Compressor.PyTree as Compressor
    psi = time*omega*180./numpy.pi
    p = phase(time, omega, dpsi)
    #if Cmpi.rank == 0: print(p)
    if ts is None: 
        if not os.access("tstat_%d.cgns"%p,os.F_OK):
            ts = FastS.createStatNodes(t, vars=['cylindrique'])
            Internal._createUniqueChild(ts, 'phase', 'DataArray_t', value=p)
            _compactStats(ts)
        else:
            #ts=Fast.loadFile("tstat_%d.cgns"%p,mpirun=True)
            ts = Cmpi.convertFile2PyTree("tstat_%d.cgns"%p, proc=Cmpi.rank)
            _compactStats(ts)
    else:
        phase_ts = Internal.getNodeFromName1(ts,'phase')
        phase_ts = Internal.getValue(phase_ts)
        if phase_ts != p:
            Compressor._compressCartesian(ts)
            Cmpi.convertPyTree2File(ts, 'tstat_%d.cgns'%phase_ts)
            if os.access("tstat_%d.cgns"%p, os.F_OK):
                ts = Cmpi.convertFile2PyTree("tstat_%d.cgns"%p, proc=Cmpi.rank)
                _compactStats(ts)
            else : 
                ts = FastS.createStatNodes(t,vars=['cylindrique'])
                Internal._createUniqueChild(ts,'phase','DataArray_t', value=p)
                _compactStats(ts)

    _computeStats(t, ts, metrics)
    return ts

#==============================================================================
# compute enstropy et TKE in place
#==============================================================================
def _computeEnstrophy(t, metrics, time):
    own   = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy


    zones = Internal.getZones(t)
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)

    (enstrophie, tke) = fasts.computePT_enstrophy(zones,  metrics, FastC.HOOK )

    return (enstrophie, tke)

#==============================================================================
# Post: compute variable in place
#==============================================================================
def _computeVariables(t, metrics, varlist, order=2):
    """Compute given variables."""


    if   isinstance(varlist, str): vars = [varlist]
    elif isinstance(varlist, list): vars = varlist
    else: raise ValueError("_computeVariables: last argument must be a string or list of strings.")

    lcompact_Q    = False
    lcompact_Enst = False
    lcompact_drodt= False
    lcompact_Rotx = False

    own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    flag = 0
    var_loc = []

    for var in vars:
       #print 'var=',var
       if var == 'QCriterion':
          flag += 1
          var_loc.append('QCriterion')
       elif var == 'QpCriterion':
          flag += 2
          var_loc.append('QpCriterion')
       elif var == 'Enstrophy':
          flag += 10
          var_loc.append('Enstrophy')
       elif var == 'RotX':
          flag += 100
          var_loc.append('RotX')
       elif var == 'dDensitydt':
          flag += 1000
          var_loc.append('dDensitydt')

    #print flag,var_loc
    ##verifie si les noeux existe dans l'arbre
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        size = (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        solution = Internal.getNodeFromName1(z, 'FlowSolution#Centers')

        for var in var_loc:
            node = Internal.getNodeFromName1(solution, var)
            if node is None:
                 if var=='QCriterion' or var=='QpCriterion': lcompact_Q     = True
                 if var=='RotX'                            : lcompact_Rotx  = True
                 if var=='Enstrophy'                       : lcompact_Enst  = True
                 if var=='dDensitydt'                      : lcompact_drodt = True
                 tmp = numpy.ones( (dim[1]-1,dim[2]-1,dim[3]-1) , dtype=numpy.float64)
                 Internal.createChild(solution, var, 'DataArray_t', tmp)

    if var_loc != []:
       if lcompact_Q:     FastC._compact(zones, fields=['centers:QCriterion'], dtloc=dtloc)
       if lcompact_Enst:  FastC._compact(zones, fields=['centers:Enstrophy' ], dtloc=dtloc)
       if lcompact_Rotx:  FastC._compact(zones, fields=['centers:RotX'      ], dtloc=dtloc)
       if lcompact_drodt: FastC._compact(zones, fields=['centers:dDensitydt'], dtloc=dtloc)

       own   = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
       dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
       dtloc = Internal.getValue(dtloc)                       # tab numpy

       # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
       if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)
       fasts.computePT_variables(zones,  metrics, FastC.HOOK, flag, order)
    return None

#==============================================================================
# Post: compute gradient in place
#==============================================================================
def _computeGrad(t, metrics, varlist, order=2):
    """Compute gradient of fiven variables."""


    if isinstance(varlist, str): vars = [varlist]
    elif isinstance(varlist, list): vars = varlist
    else: raise ValueError("_computeGrad: last argument must be a string or list of strings.")

    own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    var_zones = []
    vargrad   = []
    lgrad     = False

    ##verifie si les noeuds existent dans l'arbre
    zones = Internal.getZones(t)
    for z in zones:
        #
        dim   = Internal.getZoneDim(z)
        size = (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        solution = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
        cgrad=0
        var_loc = []
        for var in vars:
            node = Internal.getNodeFromName1(solution, var)
            #verifie si la variable existe dans l'arbre
            if node is None:
               print("no", var, "in tree: gradient is not computed.")
            else:
               var_loc.append(var)
               lgrad    = True
               lcompact = False
               vargrad.append('gradx'+var)
               vargrad.append('grady'+var)
               vargrad.append('gradz'+var)
               #creation noeuds gradient si necessaire
               for grad in ['gradx'+var, 'grady'+var, 'gradz'+var]:
                   node_grad = Internal.getNodeFromName1(solution, grad)
                   if node_grad is None:
                      #tmp = numpy.empty( (dim[1]-1,dim[2]-1,dim[3]-1) , dtype=numpy.float64)
                      tmp = numpy.ones( (dim[1]-1,dim[2]-1,dim[3]-1) , dtype=numpy.float64)
                      Internal.createChild(solution, grad, 'DataArray_t', tmp)
                      lcompact = True

               if lcompact: FastC._compact(z, fields=['centers:'+vargrad[cgrad],'centers:'+vargrad[cgrad+1],'centers:'+vargrad[cgrad+2]], dtloc=dtloc)
               cgrad+=3

        var_zones.append(var_loc)
      
    if lgrad:
       dtloc = Internal.getNodeFromName2(t , '.Solver#dtloc')  # noeud
       dtloc = Internal.getValue(dtloc)                        # tab numpy

       # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
       if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)
       fasts.computePT_gradient(zones,  metrics, var_zones, vargrad, FastC.HOOK, order)
    return None

#==============================================================================
# Display
#==============================================================================
def displayTemporalCriteria(t, metrics, nitrun, format=None, gmres=None, verbose='firstlast', stopAtNan=True):
    """Display CFL and convergence information."""
    return display_temporal_criteria(t, metrics, nitrun, format, gmres, verbose, stopAtNan)

def display_temporal_criteria(t, metrics, nitrun, format=None, gmres=None, verbose='firstlast', stopAtNan=True):
    own          = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc        = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
    dtloc_numpy  = Internal.getValue(dtloc)
    nssiter      = int(dtloc_numpy[0])

    zones    = Internal.getZones(t)
    nzones	 = len(zones)

    #a = Internal.getNodeFromName2(zones[0], 'model')
    #model = Internal.getValue(a)
    #neq = 5
    #if (model == 'nsspalart' or model =='NSTurbulent'): neq = 6

    neq_max = 6

    cvg_numpy = numpy.empty((nzones,4*neq_max), dtype=numpy.float64)
    # sortie sur stdout "simple precision"
    lft = 1
    # sortie sur stdout "double precision"
    if format == "double": lft = 0
    # sortie sur Fichier Fortran binaire
    elif format == "flush": lft = -1
    # sortie sur Fichier Fortran binaire
    elif format == "ascii": lft = -2
    # enregistrement dans l'arbre uniquement
    elif format == "store": lft = 3

    iverb = 0
    if verbose != 'firstlast': iverb=2
    residu = fasts.display_ss_iteration(zones, metrics, cvg_numpy, nitrun, nssiter, lft, iverb)

    if stopAtNan:
        finite = Cmpi.isFinite(t, var='centers:Density')
        if not finite:
            import sys; sys.exit(1) # return 1 to shell

    if gmres is None: return None
    else: return residu

#==============================================================================
# IN: d: container
# IN: keys: les cles possibles
#==============================================================================
def checkKeys(d, keys):
    for i in d[2]:
        if i[0] not in keys:
            print('Warning: FastS: keyword %s is invalid.'%i[0])

#==============================================================================
# init velocity (ALE)
#==============================================================================
def _computeVelocityAle(t, metrics):

    own   = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy

    zones = Internal.getZones(t)
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it

    node = Internal.getNodeFromName1(t, '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = FastC.OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    fasts.computePT_velocity_ale(zones,  metrics, FastC.HOOK, ompmode)
    return None

#==============================================================================
# init velocity (ALE) from DADS maillage deformable
#==============================================================================
def copy_velocity_ale(t, metrics, it=0):

    own   = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy

    zones = Internal.getZones(t)
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it

    node = Internal.getNodeFromName1(t   , '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = FastC.OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    #on filtre les zones non deformables
    zones_aleDeformation=[]; metric_aleDeformation=[]
    c=0
    for z in zones:
      tmp = Internal.getNodeFromName1(z, 'Motion')
      if tmp is not None:
          zones_aleDeformation.append(z)
          metric_aleDeformation.append(metrics[c])
      c += 1

    fasts.copy_velocity_ale(zones_aleDeformation, metric_aleDeformation, FastC.HOOK, ompmode, it)
    return None

#==============================================================================
# move mesh (ALE)
#==============================================================================
def _movegrid(t):
    zones = Internal.getZones(t)
    fasts._movegrid(zones)
    return None

#==============================================================================
# Convergence
#==============================================================================
def createConvergenceHistory(t, nrec):
    """Create a node in tree to store convergence history."""
    varsR   = ['RSD_L2','RSD_oo','RSD_L2_diff','RSD_oo_diff']
    bases   = Internal.getNodesFromType1(t, 'CGNSBase_t')
    curIt   = 0
    for b in bases:
       Internal.createUniqueChild(b, 'GlobalConvergenceHistory',
                                  'ConvergenceHistory_t', value=curIt)
    zones = Internal.getZones(t)
    neq = 5
    a = Internal.getNodeFromName3(t, 'GoverningEquations')
    if a is not None:
       model = Internal.getValue(a)
       if model == 'nsspalart' or model =='NSTurbulent': neq = 6
    for z in zones:
        c = Internal.createUniqueChild(z, 'ZoneConvergenceHistory',
                                       'ConvergenceHistory_t', value=curIt)
        tmp = numpy.zeros((nrec), numpy.int32)
        Internal.createChild(c, 'IterationNumber', 'DataArray_t', tmp)
        for var in varsR:
            tmp = numpy.zeros((nrec*neq), numpy.float64)
            Internal.createChild(c, var ,'DataArray_t', tmp)
    return None

#==============================================================================
# extraction des residus du Pytree dans fichier tecplot ascii
# ==============================================================================
def write_plt_format(t, i, FileCvg, nd, it=[], RSD_L2=[], RSD_oo=[], RSD_L2_diff=[], RSD_oo_diff=[], 
                     a=[], convergence_name='ZoneConvergenceHistory'):
    node    = Internal.getNodeFromPath(t, i+'/'+convergence_name)
    lastRec = node[1][0]

    add_local = '/'+convergence_name+'/'
    it = Internal.getNodeFromPath(t,i+add_local+"IterationNumber")
    
    a      = Internal.getNodeFromPath(t, i+add_local+"RSD_L2")
     
    nrec = it[1].size
    neq  = a[1].size//nrec

    RSD_L2 = numpy.reshape(a[1],(neq,nrec),order='F')
    
    a = Internal.getNodeFromPath(t, i+add_local+"RSD_oo")
    RSD_oo = numpy.reshape(a[1],(neq,nrec),order='F')
    
    a = Internal.getNodeFromPath(t, i+add_local+"RSD_L2_diff")
    RSD_L2_diff = numpy.reshape(a[1],(neq,nrec),order='F')
    
    a = Internal.getNodeFromPath(t, i+add_local+"RSD_oo_diff")
    RSD_oo_diff = numpy.reshape(a[1],(neq,nrec),order='F')
    
    a='"'
    if neq == 5:
        var="VARIABLES = it RO_l2 ROU_l2 ROV_l2 ROW_l2 ROE_l2 ROoo ROUoo ROVoo ROWoo ROEoo RO_l2_diff ROU_l2_diff ROV_l2_diff ROW_l2_diff ROE_l2_diff ROoo_diff ROUoo_diff ROVoo_diff ROWoo_diff ROEoo_diff\n"
    if neq == 6:
        var="VARIABLES = it RO_l2 ROU_l2 ROV_l2 ROW_l2 ROE_l2 NUT_l2 ROoo ROUoo ROVoo ROWoo ROEoo NUToo RO_l2_diff ROU_l2_diff ROV_l2_diff ROW_l2_diff ROE_l2_diff NUT_l2_diff ROoo_diff ROUoo_diff ROVoo_diff ROWoo_diff ROEoo_diff NUToo_diff\n"
    if nd == 0: FileCvg.write("%s"%(var))
    nd =+ 1
    FileCvg.write("ZONE T=%sbloc %s %s I=%d F=POINT\n"%(a,i,a,lastRec))
    c = it[1]
    for l in range(lastRec):
       a  = ""
       for k in range(neq): a = a+"{0:7f} ".format(RSD_L2     [(k,l)])
       for k in range(neq): a = a+"{0:7f} ".format(RSD_oo     [(k,l)])
       for k in range(neq): a = a+"{0:7f} ".format(RSD_L2_diff[(k,l)])
       for k in range(neq): a = a+"{0:7f} ".format(RSD_oo_diff[(k,l)])           
       FileCvg.write('%s %s\n'%(1+c[l],a))
    return nd


# IN: t: tree with ConvergenceHistory ()
# IN: fileout: fileName for output of residuas in tp format
# IN: perZones: write ConvergenceHistory for each zones
# IN: perBases: write ConvergenceHistory for each base
def _extractConvergenceHistory(t, fileout, perZones=True, perBases=True):
    """Extract residuals in an ascii file (tp format)."""
    zones = Internal.getZonePaths(t)
    nd = 0
    fileCvgName = fileout
    fileCvg = open(fileCvgName, 'w')

    # Ecrit le residu de chaque zone    
    if perZones:
        for z in zones:
            nd = write_plt_format(t, z, fileCvg, nd, it=[], RSD_L2=[], RSD_oo=[], RSD_L2_diff=[], RSD_oo_diff=[], 
                                  a=[], convergence_name='ZoneConvergenceHistory')
    
    # Ecrit le residu par base
    if perBases:
        for i in Internal.getPathsFromType(t, 'CGNSBase_t'):
            base = Internal.getNodeFromPath(t, i)
            zone_check = Internal.getNodeByName(base, 'GlobalConvergenceHistory')
            if Internal.getChildren(zone_check): # if list is not empty
                nd = write_plt_format(t, i, fileCvg, nd, it=[], RSD_L2=[], RSD_oo=[], RSD_L2_diff=[], RSD_oo_diff=[], 
                                      a=[], convergence_name='GlobalConvergenceHistory')
    
    fileCvg.close()

#==============================================================================
# Cree un arbre Stress pour calcul effort
# IN: t: tree
# IN: BC: list des BC concernees ['BCWall', BCFarfield',...]
# IN: windows: ou range de window
# OUT: return arbre stress
#==============================================================================
def createStressNodes(t, BC=None, windows=None):
    """Create nodes to store stress data."""
    import Converter.GhostCells as Ghost
    try: import Transform.PyTree as T
    except: raise ImportError("createStressNodes: requires transform module.")

    vars0 = ['CoordinateX','CoordinateY','CoordinateZ','centers:cellN']
    var   = ['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity']

    bases = Internal.getBases(t)
    teff = C.newPyTree([b[0] for b in bases])

    familyBCDict = C.getFamilyBCNamesDict(t)

    own   = Internal.getNodeFromName1(t,   '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')    # noeud

    teff[2].append(own)

    zones = []
    no_z = 0
    for b0 in Internal.getNodesFromType1(t, 'CGNSBase_t'):

        b = Internal.getNodeFromName1(teff, b0[0])

        dimbase = b0[1][0]
           
        count = -1
        zones = Internal.getNodesFromType1(b0, 'Zone_t')

        for z in b0[2]:
        
            count += 1
            if z[3] != 'Zone_t': continue

            # gestion zone2d
            dimzone = Internal.getZoneDim(z)
            inck0 = 2 ; inck1 = 1
            #vargrad=['gradxVelocityX','gradyVelocityX','gradzVelocityX','gradxVelocityY','gradyVelocityY','gradzVelocityY','gradxVelocityZ','gradyVelocityZ','gradzVelocityZ','gradxTemperature','gradyTemperature','gradzTemperature', 'CoefPressure','ViscosityMolecular']
            vargrad=['gradxVelocityX','gradyVelocityX','gradzVelocityX','gradxVelocityY','gradyVelocityY','gradzVelocityY','gradxVelocityZ','gradyVelocityZ','gradzVelocityZ','gradxTemperature','gradyTemperature','gradzTemperature', 'CoefPressure','ViscosityMolecular','Density2','Pressure']
            if dimzone[3] == 2: 
                inck0=0 ; inck1=0
                #vargrad=['gradxVelocityX','gradyVelocityX','gradxVelocityY','gradyVelocityY','gradxTemperature','gradyTemperature','CoefPressure','ViscosityMolecular']
                vargrad=['gradxVelocityX','gradyVelocityX','gradxVelocityY','gradyVelocityY','gradxTemperature','gradyTemperature','CoefPressure','ViscosityMolecular','Density2','Pressure']

            varc  = []
            for i in var: varc.append('centers:'+i)
            for i in vargrad: varc.append('centers:'+i)

            if BC is not None:
                bc = Internal.getNodeFromName1(z, 'ZoneBC')
                list_bc = []
                if bc is not None: list_bc = bc[2]
            else:
                list_bc =['window']
            ific = 2
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')
            if param_int is not None: ific = param_int[1][3]
            ndf = 0
            for v in list_bc:

                selectBC = False 
                if BC is not None:
                    # nom de la BC
                    name = Internal.getValue(v)
                    if name == 'FamilySpecified': # name is a familyName
                        familyName = Internal.getNodeFromName1(v, 'FamilyName')
                        if familyName is not None:
                            familyName = Internal.getValue(familyName)
                            if familyName in BC: selectBC = True
                            if familyBCDict[familyName] in BC: selectBC = True
                    else: # name is a type
                        if name in BC: selectBC = True
                    #print('selection', name, selectBC)
                else: name = v[0] 

                if windows is None: windows = [None]

                inum = 0

                for window in windows:

                    if (selectBC) or (window is not None and z[0]==window[0]):
                        if BC is not None:
                            ptrg = Internal.getNodeFromName(v, "PointRange")
                        else: 
                            wrange = numpy.zeros(6, numpy.int32 ).reshape(3,2)
                            wrange[0,0]= window[1]
                            wrange[1,0]= window[3]
                            wrange[2,0]= window[5]
                            wrange[0,1]= window[2]
                            wrange[1,1]= window[4]
                            wrange[2,1]= window[6]
                            ptrg = ['PointRange', wrange, [], 'IndexRange_t']
                    
                        dim = Internal.getValue(ptrg)
                    
                        #dir=0,...5
                        idir = Ghost.getDirection__(dimbase, [ptrg])
                        if windows[0] is not None:
                            if   window[7] == 'imin': idir =0
                            elif window[7] == 'imax': idir =1
                            elif window[7] == 'jmin': idir =2
                            elif window[7] == 'jmax': idir =3
                            elif window[7] == 'kmin': idir =4
                            elif window[7] == 'kmax': idir =5
                    
                        inci =0 ; incj =0 ; inck =0
                        if BC is not None:
                            if   idir==0: inci= ific
                            elif idir==1: inci=-ific
                            elif idir==2: incj= ific
                            elif idir==3: incj=-ific
                            elif idir==4: inck= inck0
                            elif idir==5: inck=-inck0

                        ni = dim[0,1]-dim[0,0]
                        nj = dim[1,1]-dim[1,0]
                        nk = dim[2,1]-dim[2,0]

                        ideb = dim[0,0]+inci
                        ifin = dim[0,1]+inci
                        jdeb = dim[1,0]+incj
                        jfin = dim[1,1]+incj
                        kdeb = dim[2,0]+inck
                        kfin = dim[2,1]+inck

                        #print('subZ', ideb,ifin,jdeb,jfin,kdeb,kfin)
                        zp = T.subzone(z, (ideb,jdeb,kdeb), (ifin,jfin,kfin))
                        zp[0] = z[0]+'_'+v[0]+str(inum)
                        inum  += 1
                        C._extractVars(zp, vars0)
                        c = 0
                        for v1 in varc: 
                            C._initVars(zp, v1, c)
                            c += 1
                        #print('zone name=',zp[0])
                        #print('dim0=',dim[0,0],dim[0,1],dim[1,0],dim[1,1],dim[2,0],dim[2,1])

                        # Recuperation du GridCoordinates#Init si il existe
                        ci = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
                        if ci is not None:
                            cin = Internal.createUniqueChild(zp, 'GridCoordinates#Init', 'GridCoordinates_t')
                            cin[2] = []
                            n = Internal.getNodeFromName1(ci, 'CoordinateX')
                            sl = n[1][ideb-1:ifin, jdeb-1:jfin, kdeb-1:kfin]
                            Internal.newDataArray('CoordinateX', value=sl, parent=cin)
                            n = Internal.getNodeFromName1(ci, 'CoordinateY')
                            sl = n[1][ideb-1:ifin, jdeb-1:jfin, kdeb-1:kfin]
                            Internal.newDataArray('CoordinateY', value=sl, parent=cin)
                            n = Internal.getNodeFromName1(ci, 'CoordinateZ')
                            sl = n[1][ideb-1:ifin, jdeb-1:jfin, kdeb-1:kfin]
                            Internal.newDataArray('CoordinateZ', value=sl, parent=cin)

                        if idir <= 1:
                            ni1  = nj
                            nj1  = nk
                            jmin = jdeb -ific      #jmin
                            jmax = jfin -ific-1    #jmax
                            kmin = kdeb -inck0     #kmin
                            kmax = kfin -inck0-1   #kmax
                            imin = ideb -ific
                            imax = imin
                        elif idir <= 3:   
                            ni1  = ni
                            nj1  = nk
                            imin = ideb -ific      #imin
                            imax = ifin -ific-1    #imax
                            kmin = kdeb -inck0     #kmin
                            kmax = kfin -inck0-1   #kmax
                            jmin = jdeb -ific
                            jmax = jmin
                        else:
                            ni1  = ni
                            nj1  = nj
                            imin = ideb -ific      #imin
                            imax = ifin -ific-1    #imax
                            jmin = jdeb -ific      #jmin
                            jmax = jfin -ific-1    #jmax
                            kmin = kdeb -inck0
                            kmax = kmin

                        param_int = Internal.getNodeFromName2(zp, 'Parameter_int')  # noeud
                        if param_int is None:
                            raise ValueError("_createStressNodes: Parameter_int is missing for zone %s."%zp[0])

                        #print('loop effort',imin,imax,jmin,jmax,kmin,kmax)
                        #modifie le moeud param_int de l'arbre teffort (issue de subzone) pour la fonction initNuma
                        param_int[1]    = numpy.copy(param_int[1])    # sinon tableau partager avec param de la zone NS
                        #NIJK  pour initNuma
                        param_int[1][ 0] = ni1
                        param_int[1][ 1] = nj1
                        param_int[1][ 2] = 1
                        param_int[1][ 3] = ific
                        param_int[1][ 4] = 0
                        #IJKV  pour initNuma
                        param_int[1][20] = max(ni1-2*ific, 1)
                        param_int[1][21] = max(nj1-2*ific, 1)
                        param_int[1][22] = 1
                        param_int[1][36] = len(varc)             #NEQ    
                        param_int[1][41] = ni1*nj1               #NDIMDX
                        param_int[1][66] = 0                     #SHIFTVAR
                        #fenetre calcul flux dans arbre NS
                        param_int[1][23] = imin
                        param_int[1][24] = imax
                        param_int[1][25] = jmin
                        param_int[1][26] = jmax
                        param_int[1][27] = kmin
                        param_int[1][28] = kmax
                        #adresse stockage flu:  adr = 1 + (i-i0) +(j-j0)*Ni0 + (k-k0)*Ni0*Nj0
                        param_int[1][29] = imin                             #i0
                        param_int[1][30] = jmin                             #j0
                        param_int[1][31] = kmin                             #k0
                        param_int[1][32] = imax-imin+1                      #Ni0
                        param_int[1][33] = param_int[1][32]*(jmax-jmin+1)   #Ni0*Nj0
                        #no de la zone pour recuperer pointer zone
                        param_int[1][34] = no_z               
                        param_int[1][35] = idir+1
                        #param_int[1][36] = ndf
                    
                        #print('dim1=',imin,imax,jmin,jmax,kmin,kmax)
                        #print('idir=',idir+1)

                        FastC._compact(zp, fields=varc, init=False, dtloc=dtloc) #allocation compact uniquememnt; init dans er calcul effort

                        b[2].append(zp)

            no_z +=1
            ndf  +=1

        # ajoute le refState (if any)
        ref = Internal.getNodeFromName1(b0, 'ReferenceState')
        if ref is not None: b[2].append(ref)
        ref = Internal.getNodeFromName1(b0, 'FlowEquationSet')
        if ref is not None: b[2].append(ref)

    Internal._rmNodesByType(teff, 'ZoneGridConnectivity_t')
    Internal._rmNodesByType(teff, 'ZoneBC_t')
    Internal._rmNodesByType(teff, 'Rind_t')
    Internal._rmNodesByName(teff, '.Solver#define')
    Internal._rmNodesByName(teff, 'Parameter_real')
    Internal._rmNodesByName(teff, 'CFL_minmaxmoy')
    Internal._rmNodesByName(teff, 'type_zone')
    Internal._rmNodesByName(teff, 'model')
    Internal._rmNodesByName(teff, 'temporal_scheme')
    Internal._rmNodesByName(teff, 'FlowSolution')
    Internal._rmNodesByName(teff, 'Motion')
    
    #Internal._rmNodesByName(teff, 'GridCoordinates#Init')
    
    return teff

#==============================================================================
#
# Calcul des efforts (in place)
#
#==============================================================================
def _computeStress(t, teff, metrics, xyz_ref=(0.,0.,0.)):
    """Compute efforts in teff."""

    zones     = Internal.getZones(t)
    zones_eff = Internal.getZones(teff)
    
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    if FastC.HOOK is None:
            own    = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
            dtloc  = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
            dtloc  = Internal.getValue(dtloc)                       # tab numpy
            FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)
            nitrun =0; nstep =1

            hook1 = FastC.HOOK.copy()
            hook1.update(FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, 0) )

    else:  hook1  = FastC.HOOK

    effort  = numpy.zeros(11, numpy.float64)
    pos_eff = numpy.empty( 3, numpy.float64)
    pos_eff[0] = xyz_ref[0]; pos_eff[1] = xyz_ref[1]; pos_eff[2] = xyz_ref[2]

    fasts.compute_effort(zones, zones_eff, metrics, hook1, effort, pos_eff)

    return effort

#==============================================================================
#
# IBC ordre 0: preparation
# surf: arbre avec une base par obstacle (plusieurs zones possibles)
#==============================================================================
def setIBCData_zero(t, surf, dim=None):

  bases = Internal.getBases(t)
  zones = Internal.getZones(t)

  # recup des obstacles
  bodies = []
  for b in Internal.getBases(surf):
      body = Internal.getZones(b)
      bodies.append(body)

  # Matrice de masquage (arbre d'assemblage)
  nbodies = len(bodies)
  nbases = len(bases)
  BM = numpy.ones((nbases,nbodies), dtype=numpy.int32)
  if dim is None:
     sol= Internal.getNodeFromName1(zones[0],'FlowSolution#Centers')
     nk = Internal.getNodeFromName1(sol,'Density')[1].shape[2]
     if nk ==1 : dim = 2
     else:       dim = 3

  C._initVars(t,"{centers:cellN_IBC}=1.")
  #X._blankCells(t, bodies, BM, blankingType="center_in", dim=dim, delta=0., cellNName='cellN_IBC')
  t = X.blankCells(t, bodies, BM, blankingType="center_in", dim= dim,delta=0., cellNName='cellN_IBC')

  numz = {}
  zones = Internal.getZones(t)
  for z in zones:
     cellN_node = Internal.getNodeFromName2(z,'cellN_IBC')
     #cellN_node[0]= 'cellN_IBC'
     cellN_IBC = cellN_node[1]

     Internal._rmNodesByName(z, "IBC")

     if 0 in cellN_IBC:
        ibc = numpy.ones( 7, dtype=numpy.int32)
        itemindex = numpy.where( cellN_IBC==0 )
        ibc[0]=1
        #formule valable pour 2Ghostcells: sinon -ific+1
        ibc[1]= numpy.amin( itemindex[0] )-1
        ibc[2]= numpy.amax( itemindex[0] )-1
        ibc[3]= numpy.amin( itemindex[1] )-1
        ibc[4]= numpy.amax( itemindex[1] )-1

        if cellN_IBC.shape[2]!=1:
           #print'verifier en 3d'
           ibc[5]= numpy.amin( itemindex[2] )-1
           ibc[6]= numpy.amax( itemindex[2] )-1

        #on retranche 1 pour avoir une bonne visu et un terme src plus compact
        cellN_IBC[0:] -= 1

        print('Zone', z[0],': plage IBC=', ibc[1:7])

        numz["IBC"]  = ibc

        FastC._setNum2Zones(z, numz)
     else:
        C._rmVars(z, ['centers:cellN_IBC'])

  return t

#==============================================================================
#
# rmGhost robuste 2d et 3d
#==============================================================================
def rmGhostCells(t, depth, adaptBCs=1):

 try: import Transform.PyTree as T
 except: raise ImportError("rmGhostCells: requires Transform module.")

 zones = Internal.getZones(t)

 basename = Internal.getNodeFromType1(t, 'CGNSBase_t')[0]

 zout=[]
 for z in zones:
   dim = Internal.getZoneDim(z)
   zname = z[0]
   if dim[3]==2:  #zone 2D
      coordz = Internal.getNodeFromName2( z, "CoordinateZ" )[1]
      dz = coordz[2,2,1]-coordz[2,2,0]
      z = T.subzone(z,(1,1,1),(-1,-1,1))
      z = Internal.rmGhostCells(z, z, 2, adaptBCs=adaptBCs)
      z = T.addkplane(z,1)

      coordz = Internal.getNodeFromName2( z, "CoordinateZ" )[1]
      b =z[1]
      jmax = b[1,1]+1
      imax = b[0,1]+1
      coordz[0:imax,0:jmax,1] = coordz[2,2,0] +dz
   else:
      z = Internal.rmGhostCells(z, z, 2, adaptBCs=adaptBCs)

   z[0] = zname
   zout.append(z)

 tp = C.newPyTree([basename])
 # Attache les zones
 tp[2][1][2] += zout

 return tp

#==============================================================================
# display CPU efficiency diagnostic
#==============================================================================
def display_cpu_efficiency(t, mask_cpu=0.08, mask_cell=0.01, diag='compact', FILEOUT='listZonesSlow.dat', FILEOUT1='diagCPU.dat', RECORD=None):

 own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
 dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

 node = Internal.getNodeFromName1(t, '.Solver#define')
 node = Internal.getNodeFromName1(node, 'omp_mode')
 ompmode = FastC.OMP_MODE
 if  node is not None: ompmode = Internal.getValue(node)

 dtloc = Internal.getValue(dtloc) # tab numpy
 ss_iteration  = int(dtloc[0])

 timer_omp = FastC.HOOK["TIMER_OMP"]
 ADR = OMP_NUM_THREADS*2*(ss_iteration)
 echant    =  timer_omp[ ADR ]
 if echant == 0.:
   print('nombre iterations insuffisant pour diagnostic: nitrun * ss_iteration > 15')
   return None

 cellCount = numpy.zeros(2*OMP_NUM_THREADS, dtype=numpy.int32)

 zones = Internal.getZones(t)
 tps_percell     =0.
 tps_percell_max =0.
 tps_percell_min =0.
 cout            =0.
 cells_tot       =0
 data            ={}
 f1 = open(FILEOUT1, 'w')
 for z in zones:
    echant    =  timer_omp[ ADR ]
    param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
    solver_def= Internal.getNodeFromName2(z, '.Solver#define')
    if ompmode == 1:
       Ptomp       = param_int[69]
       PtrIterOmp  = param_int[Ptomp]
       PtZoneomp   = param_int[PtrIterOmp]
       NbreThreads = param_int[ PtZoneomp  + OMP_NUM_THREADS ]
    else:
       NbreThreads = OMP_NUM_THREADS

    if diag == 'compact':
      tps_zone_percell     =0.
      tps_zone_percell_max =0.
      tps_zone_percell_min =100000.
      ithread_max          = 0
      ijkv                 = param_int[20]*param_int[21]*param_int[22]
      for i in range(OMP_NUM_THREADS):
          if param_int[34]==0:#non ALE
              cellCount[2*i]+=timer_omp[ ADR + 2+i*2 ]
          else:
              cellCount[2*i+1]+=timer_omp[ ADR + 2+i*2 ]

          if ompmode == 1:
            ithread = param_int[ PtZoneomp  +  i ]
            if ithread != -2:
               #print "check", z[0],timer_omp[ ADR + 1+i*2 ]/echant, timer_omp[ ADR + 2+i*2 ], i,"echant=",echant,ijkv
               #tps_zone_percell += timer_omp[ ADR + 1+ithread ]
               tps_zone_percell += timer_omp[ ADR + 1+i*2 ]*timer_omp[ ADR + 2+i*2 ]/float(ijkv)*NbreThreads   #tps * Nb cell
               if timer_omp[ ADR + 1+i*2 ] > tps_zone_percell_max: 
                    tps_zone_percell_max = timer_omp[ ADR + 1+i*2 ]
                    ithread_max          = ithread
               if timer_omp[ ADR + 1+i*2 ] < tps_zone_percell_min: 
                    tps_zone_percell_min = timer_omp[ ADR + 1+i*2 ]
                    ithread_min          = ithread
          else:
            tps_zone_percell += timer_omp[ ADR + 1+i*2 ]*timer_omp[ ADR + 2+i*2 ]/float(ijkv)*NbreThreads
            if timer_omp[ ADR + 1+i*2 ] > tps_zone_percell_max: 
                 tps_zone_percell_max = timer_omp[ ADR + 1+i*2 ]
                 ithread_max          = i
            if timer_omp[ ADR + 1+i*2 ] < tps_zone_percell_min: 
                 tps_zone_percell_min = timer_omp[ ADR + 1+i*2 ]
                 ithread_min          = i
      
      tps_percell     += tps_zone_percell/echant/NbreThreads*ijkv
      tps_percell_max += tps_zone_percell_max/echant*ijkv
      tps_percell_min += tps_zone_percell_min/echant*ijkv
      cout_zone        = tps_zone_percell/echant/NbreThreads*ijkv
      #tps_percell     += tps_zone_percell/echant/NbreThreads/
      #tps_percell_max += tps_zone_percell_max/echant
      #cout_zone        = tps_zone_percell/echant/NbreThreads
      cout            += cout_zone
      data[ z[0] ]   = [ tps_zone_percell/echant/NbreThreads, cout_zone]

      f1.write('cpumoy/cell=  '+str(tps_zone_percell/echant/NbreThreads)+"  cpumaxmin= "+str(tps_zone_percell_max/echant)+str(tps_zone_percell_min/echant) +" th lent/rap= "+str(ithread_max)+str(ithread_min)+"  Nthtread actif="+str(NbreThreads)+" dim zone="+str(ijkv)+" dim tot="+str(cells_tot)+"  "+z[0]+'\n')

      cells_tot   += ijkv

      if RECORD is None: print('cpu/cell zone=', z[0],'moy=',tps_zone_percell/echant/NbreThreads,'maxmin=',tps_zone_percell_max/echant, tps_zone_percell_min/echant,'th maxmin=', ithread_max, ithread_min, 'Nbtread actif=', NbreThreads, ' dim zone=',ijkv,'dim tot=',cells_tot)
      if RECORD is not None:
          tape = tps_zone_percell/echant/NbreThreads

      tps_zone_percell = max(tps_zone_percell, 1.e-11)
      tps_zone_percell_max = max(tps_zone_percell_max, 1.e-11)

      perfo = numpy.empty(2, dtype=numpy.float64)
      perfo[0]= int(echant*NbreThreads/tps_zone_percell)
      perfo[1]= int(echant/tps_zone_percell_max)
      ## CUPS  : [0] MEAN VALUE [1] MIN VALUE
      ## 1/CUPS=mus/iter(time * sub-iter)/cellspercore

      Internal.createUniqueChild(solver_def, 'Cups', 'DataArray_t', value=perfo)

    else:
      for i in range(OMP_NUM_THREADS):
          if ompmode ==1:
            ithread = param_int[ PtZoneomp  +  i ]
            if ithread != -2:
               if RECORD is None: print('zone= ', z[0], 'cpu= ',timer_omp[ ADR + 1+ithread*2 ]/echant,' th=  ', ithread, 'echant= ', echant)
          else:
            ithread = i
            if RECORD is None: print('zone= ', z[0], 'cpu= ',timer_omp[ ADR + 1+ithread*2 ]/echant,' th=  ', ithread, 'echant= ', echant)

    ADR+= OMP_NUM_THREADS*2+1

 for i in range(OMP_NUM_THREADS):
   print('Nbr cellule Fixe/ALE: ',cellCount[2*i],cellCount[2*i+1] , "pour th=", i)

 tps_percell/=cells_tot
 if RECORD is None: print('cpu moyen %cell en microsec: ', tps_percell, tps_percell_max/cells_tot,tps_percell_min/cells_tot )

 f1.write('cpu moyen %cell en microsec: '+str(tps_percell)+'   '+str(tps_percell_max/cells_tot)+str(tps_percell_min/cells_tot)+ '\n')
 f1.close()

 f = open(FILEOUT, 'w')
 sizeZones=[]
 for z in zones:
    param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
    cout_relatif   = data[z[0]][0]/tps_percell-1.
    effort_relatif = data[z[0]][1]/cout
    if cout_relatif > mask_cpu and effort_relatif > mask_cell:
       sizeZone=[param_int[20],param_int[21],param_int[22]]
       if sizeZone not in sizeZones:
          sizeZones.append(sizeZone)
          if RECORD is None: print('zone ', z[0],'(',param_int[20],param_int[21],param_int[22],'):  surcout cpu= ', cout_relatif ,' , temps necessaire a cette zone (%)=', effort_relatif)

          f.write(z[0]+','+str(param_int[20])+","+str(param_int[21])+","+str(param_int[22])+","+str(param_int[25])+","+str(param_int[27])+","+str(param_int[29])+","+str(param_int[33])+'\n')
 f.close()

 if RECORD is not None: return tape
 else: return None

#
#====================================================================================
# Calcul info IBM conservatif pour l'arbre t
#====================================================================================
def _ConservativeWallIbm(t, tc, CHECK=False, zNameCheck='nada'):

 dico ={}
 families=[]
 c=0
 ztg=zNameCheck

 tmp = Internal.getNodeFromName(t, 'FlowEquationSet')
 #DimPb = Internal.getNodeFromName(tmp, 'EquationDimension')[1][0]
 DimPb = 3
 for zc in Internal.getZones(tc):
   #zc = zone donneuse
   bcs = Internal.getNodesFromName(zc,"IBCD*")
   for bc in bcs:
       lprint =False
       znameR = Internal.getValue(bc)
       family = Internal.getNodeFromName(bc, 'FamilyName')
       if family is not None: FamName = Internal.getValue(family)
       else: FamName='all'
       if FamName not in families: families.append(FamName)
       if znameR == ztg: lprint =True

       name_paroi = znameR + '#' + FamName
       if name_paroi not in dico:
            dico[name_paroi]=[0,0,0,0,0,0]

       if lprint: print ("zone=",zc[0], bc[0], FamName)
       #zone receveuse
       zR      = Internal.getNodeFromName(t, znameR)
       #cellN   = Internal.getNodeFromName(zR, "cellN#Init")[1]
       cellN   = Internal.getNodeFromName(zR, "cellN")[1]
       #on filtre les zones ou il n'y a pas de calcul NS; Sinon double  prise en compte CL
       if 1 in cellN:
         size    = numpy.size(cellN)
         sh      = numpy.shape(cellN)
         ptlistD = Internal.getNodeFromName(bc, "PointListDonor")[1]
         celln  = cellN.reshape((size), order='F')
         if lprint: print ("zone=",zc[0], bc[0],"shape=", sh)
         for l in range(numpy.size(ptlistD)):

            k =  ptlistD[l]//(sh[0]*sh[1])
            j = (ptlistD[l] -k*sh[0]*sh[1])//sh[0]
            i =  ptlistD[l] -k*sh[0]*sh[1] -j*sh[0]

            if lprint:
               print('croix ',celln[ ptlistD[l] ],celln[ ptlistD[l] +1 ],celln[ ptlistD[l] -1 ],celln[ ptlistD[l] + sh[0]],celln[ ptlistD[l] -sh[0] ], l, ptlistD[l], ptlistD[l] -sh[0],"ijk=",i,j,k)
            if celln[ ptlistD[l] -1 ]==1 and  ptlistD[l]-2 >= 0 :
              if lprint: print("imax", ptlistD[l])
              dico[name_paroi][1] +=1
            if celln[ ptlistD[l] +1 ]==1 and  ptlistD[l]+2 < size :
              if lprint: print("imin", ptlistD[l])
              dico[name_paroi][0] +=1
            if celln[ ptlistD[l] - sh[0] ]==1 and  ptlistD[l]-2*sh[0] >= 0 :
              if lprint: print("jmax", ptlistD[l])
              dico[name_paroi][3] +=1
            if celln[ ptlistD[l] + sh[0] ]==1 and  ptlistD[l]+2*sh[0] < size :
              if lprint: print("jmin", ptlistD[l])
              dico[name_paroi][2] +=1
            # cas 3D
            if DimPb ==3 :
              if celln[ ptlistD[l] - sh[0]*sh[1] ]==1 and  ptlistD[l]-2*sh[0]*sh[1] >= 0 :
                if lprint: print("kmax", ptlistD[l])
                dico[name_paroi][5] +=1
              if celln[ ptlistD[l] + sh[0]*sh[1] ]==1 and  ptlistD[l]+2*sh[0]*sh[1] < size :
                if lprint: print("kmin", ptlistD[l])
                dico[name_paroi][4] +=1

 #dimensinonement numpy stockage No face
 for key in dico:
   if ztg in key: print("SIZE FLUX=", dico[key][0],dico[key][1],dico[key][2],dico[key][3],dico[key][4],dico[key][5])
   #print (key, dico[key])
   zsrname = key.split('#')
   zname   = zsrname[0]
   FamName = zsrname[1]

   zR      = Internal.getNodeFromName(t, zname)
   Internal.createUniqueChild(zR, 'Conservative_Flux', 'UserDefinedData_t')
   tmp1     = Internal.getNodeFromName(zR, 'Conservative_Flux')
   Internal.createUniqueChild(tmp1 ,FamName, 'UserDefinedData_t')
   tmp      = Internal.getNodeFromName(tmp1, FamName)

   FaceList = ['Face_imin','Face_imax', 'Face_jmin','Face_jmax','Face_kmin','Face_kmax']
   if DimPb==2: FaceList = ['Face_imin','Face_imax', 'Face_jmin','Face_jmax']

   c = 0
   for face in FaceList:
       datap = numpy.zeros( dico[key][c], numpy.int32)
       Internal.createUniqueChild(tmp, face, 'DataArray_t', datap)
       c+=1

 if CHECK:
   for z in Internal.getZones(t):
      C._initVars(t,"{centers:cellI#conservatif}=0")
      C._initVars(t,"{centers:cellJ#conservatif}=0")
      if DimPb == 3: C._initVars(t,"{centers:cellK#conservatif}=0")

 #stockage No face
 for zc in Internal.getZones(tc):
   #zc = zone donneuse
   bcs = Internal.getNodesFromName(zc,"IBCD*")
   for bc in bcs:
      znameR = Internal.getValue(bc)

      family = Internal.getNodeFromName(bc, 'FamilyName')
      if family is not None: FamName = Internal.getValue(family)
      else: FamName='all'

      lprint =False
      if znameR == ztg: lprint =True
      #if not lprint  : continue
      #if FamName != 'cuve3': continue
      if lprint: print ("zone=",zc[0], "rac=", bc[0], 'fam=', FamName)
      #zone receveuse
      zR      = Internal.getNodeFromName(t, znameR)
      #cellN   = Internal.getNodeFromName(zR, "cellN#Init")[1]
      cellN   = Internal.getNodeFromName(zR, "cellN")[1]
      if CHECK: checKi  = Internal.getNodeFromName(zR, "cellI#conservatif")[1]
      if CHECK: checKj  = Internal.getNodeFromName(zR, "cellJ#conservatif")[1]
      if CHECK and DimPb == 3: checKk  = Internal.getNodeFromName(zR, "cellK#conservatif")[1]

      tmp1    = Internal.getNodeFromName(zR, "Conservative_Flux")
      tmp     = Internal.getNodeFromName(tmp1, FamName)
      Fimin   = Internal.getNodeFromName(tmp, "Face_imin")[1]
      Fimax   = Internal.getNodeFromName(tmp, "Face_imax")[1]
      Fjmin   = Internal.getNodeFromName(tmp, "Face_jmin")[1]
      Fjmax   = Internal.getNodeFromName(tmp, "Face_jmax")[1]
      sz_imin=  numpy.size(Fimin)
      sz_imax=  numpy.size(Fimax)
      sz_jmin=  numpy.size(Fjmin)
      sz_jmax=  numpy.size(Fjmax)
      sz_kmin=  0
      sz_kmax=  0

      if DimPb == 3:
        Fkmin   = Internal.getNodeFromName(tmp, "Face_kmin")[1]
        Fkmax   = Internal.getNodeFromName(tmp, "Face_kmax")[1]
        sz_kmin=  numpy.size(Fkmin)
        sz_kmax=  numpy.size(Fkmax)

      if lprint: print("     size fluxes", sz_imin,sz_imax,sz_jmin,sz_jmax,sz_kmin,sz_kmax )

      size    = numpy.size(cellN)
      sh      = numpy.shape(cellN)
      ptlistD = Internal.getNodeFromName(bc, "PointListDonor")[1]
      celln  = cellN.reshape((size), order='F')
      if CHECK:
         check_i = checKi.reshape((size), order='F')
         check_j = checKj.reshape((size), order='F')
         if DimPb == 3: check_k = checKk.reshape((size), order='F')

      val = families.index(FamName) +1

      #on skipp les zones ou il n'y a pas de points calcules
      if 1 in cellN:
        for l in range(numpy.size(ptlistD)):

           if celln[ ptlistD[l] -1 ]==1 and  ptlistD[l]-2 >= 0 :
             c                  = Fimax[ sz_imax-1]  #astuce: on stocke le compteur dans la derniere case memoire
             Fimax[c]           = ptlistD[l]+1
             if CHECK: check_i[ ptlistD[l] ] += val*100
             if c != sz_imax-1: Fimax[ sz_imax-1] += 1
             if lprint: print("            verif imax", c, sz_imax, Fimax[ sz_imax-1], 'val=', Fimax[c], ptlistD[l], l,'celln=',  celln[Fimax[c]-1], celln[Fimax[c]-2] )

           if celln[ ptlistD[l] +1 ]==1 and ptlistD[l]+2 < size :
             c                  = Fimin[ sz_imin-1]  
             Fimin[c]           = ptlistD[l]+1
             if CHECK: check_i[ ptlistD[l] ] += val*100
             if c != sz_imin-1: Fimin[ sz_imin-1] += 1
             if lprint: print("verif imin", c, sz_imin, Fimin[ sz_imin-1], Fimin[c], ptlistD[l], l,'celln=', celln[Fimin[c]-1], celln[Fimin[c]])

           if celln[ ptlistD[l] - sh[0] ]==1 and  ptlistD[l]-2*sh[0] >= 0 :
             c                  = Fjmax[ sz_jmax-1]  
             Fjmax[c]           = ptlistD[l]+1
             if CHECK: check_j[ ptlistD[l] ] += val*10000
             if c != sz_jmax-1: Fjmax[ sz_jmax-1] += 1
             if lprint: print("verif jmax", c, sz_jmax, Fjmax[ sz_jmax-1], Fjmax[c], ptlistD[l], l,'celln=', celln[Fjmax[c]-1], celln[Fjmax[c]-1-sh[0]])

           if celln[ ptlistD[l] +sh[0] ]==1 and  ptlistD[l]+2*sh[0] < size :
             c                  = Fjmin[ sz_jmin-1] 
             Fjmin[c]           = ptlistD[l]+1
             if CHECK: check_j[ ptlistD[l] ] += val*10000
             if c != sz_jmin-1: Fjmin[ sz_jmin-1] += 1
             if lprint: print("verif jmin", c, sz_jmin, Fjmin[ sz_jmin-1], Fjmin[c], ptlistD[l], l, 'celln=', celln[Fjmin[c]-1], celln[Fjmin[c]-1+sh[0]])

           if DimPb == 3:
             if celln[ ptlistD[l] - sh[0]*sh[1] ]==1 and  ptlistD[l]-2*sh[0]*sh[1] >= 0 :
                c                  = Fkmax[ sz_kmax-1]  
                Fkmax[c]           = ptlistD[l]+1
                if CHECK: check_k[ ptlistD[l] ] += val*1000000
                if c != sz_kmax-1: Fkmax[ sz_kmax-1] += 1
                if lprint: print("verif kmax", c, sz_kmax, Fkmax[ sz_kmax-1], Fkmax[c], ptlistD[l], l,'celln=', celln[Fkmax[c]-1], celln[Fkmax[c]-1-sh[0]*sh[1]])

             if celln[ ptlistD[l] +sh[0]*sh[1] ]==1 and  ptlistD[l]+2*sh[0]*sh[1] < size :
               c                  = Fkmin[ sz_kmin-1] 
               Fkmin[c]           = ptlistD[l]+1
               if CHECK: check_k[ ptlistD[l] ] += val*1000000
               if c != sz_kmin-1: Fkmin[ sz_kmin-1] += 1
               if lprint: print("verif kmin", c, sz_kmin, Fkmin[ sz_kmin-1], Fkmin[c], ptlistD[l], l, 'celln=', celln[Fkmin[c]-1], celln[Fkmin[c]-1+sh[0]*sh[1]])

#====================================================================================
# Calcule la CFL en chaque maille et la met dans l'arbre t pour dtloc instationnaire
#====================================================================================
def computeCFL_dtlocal(t, isconv=1, isvisc=1, isSound=1):

    import math
    import Transform.PyTree as T

    t = C.initVars(t,'centers:CFL', 0.)

    dimPb = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dimPb = Internal.getValue(node)

    liste_BCPeriodiques = []

    """
    BCMatch = Internal.getNodesFromType(t, 'GridConnectivity1to1_t')
    for match in BCMatch:
        perio = Internal.getNodesFromName(match, 'GridConnectivityProperty')
        if perio != [] :
            rotCenter  = Internal.getNodeFromName(perio, 'RotationCenter')[1]
            rotAngle   = Internal.getNodeFromName(perio, 'RotationAngle')[1]
            translation = Internal.getNodeFromName(perio, 'Translation')[1]

            list1 = [abs(translation).tolist(),abs(rotAngle).tolist(),rotCenter.tolist()]

            if (list1 not in liste_BCPeriodiques) :
                 liste_BCPeriodiques.append(list1)
    """
    

    (t,tc,metrics) = warmup(t, tc=None) ### Oblige d'appeler warmup afin de construire les metrics necessaires au calcul de la CFL

    zones = Internal.getZones(t)

    dim = Internal.getZoneDim(zones[0])

    print('dimPb= ', dimPb)

    fasts.prep_cfl(zones, metrics,1,1,1, isconv, isvisc, isSound)

    #t = Internal.rmGhostCells(t, t, 2)

    C._rmVars(t, 'Temperature_M1')
    C._rmVars(t, 'Temperature_P1')
    C._rmVars(t, 'Density_M1')
    C._rmVars(t, 'Density_P1')
    C._rmVars(t, 'VelocityX_M1')
    C._rmVars(t, 'VelocityX_P1')
    C._rmVars(t, 'VelocityY_M1')
    C._rmVars(t, 'VelocityY_P1')
    C._rmVars(t, 'VelocityZ_M1')
    C._rmVars(t, 'VelocityZ_P1')

    cflmax_glob = 0.0
    cflmin_glob = 1.0e15

    zones = Internal.getZones(t)

    for z in zones:
        if dimPb == 3:
           dim = Internal.getZoneDim(z)
           zp = T.subzone(z, (3,3,3), (dim[1]-2,dim[2]-2,dim[3]-2))
           cflmax = C.getMaxValue(zp, 'centers:CFL')
           cflmin = C.getMinValue(zp, 'centers:CFL')
           print('cflmin= ', cflmin)
           if cflmax > cflmax_glob:cflmax_glob = cflmax
           if cflmin < cflmin_glob:cflmin_glob = cflmin
        else:
           dim = Internal.getZoneDim(z)
           zp = T.subzone(z, (3,3,1), (dim[1]-2,dim[2]-2,1))
           cflmax = C.getMaxValue(zp, 'centers:CFL')
           cflmin = C.getMinValue(zp, 'centers:CFL')
           print('cflmin= ', cflmin)
           if cflmax > cflmax_glob:cflmax_glob = cflmax
           if cflmin < cflmin_glob:cflmin_glob = cflmin


    dtmin = 1.0/cflmax_glob
    dtmax = 1.0/cflmin_glob

    exposant_max = math.floor(math.log2(dtmax/dtmin))

    print('exposant maximal = ', exposant_max)
    print(dtmin)
    print(dtmax)

    #print(liste_BCPeriodiques)


    return t, exposant_max


#==============================================================================
# decoupe maillage pour dtloc instationnaire
#==============================================================================
##[TODO] Decision on keeping _decoupe2
##[TODO] Further testing on _decoupe2 & _decoupe4
## NOTE: _decoupe2 destroys all BCs and raccords. It is recommended to check these two
##       before proceeding.
def _decoupe2(t, exposantMax = 2, NP=0):

    import Transform.PyTree as T

    fes = Internal.getNodeFromName2(t,'FlowEquationSet')
    rs = Internal.getNodeFromName2(t,'ReferenceState')

    dtmin = 100000.
    taille_bloc = 25
    zones = []
    zones_decoupe = []

    dicoTps = {}
    numZone = []

    liste_BCPeriodiques = []

    BCMatch = Internal.getNodesFromType(t, 'GridConnectivity1to1_t')
    for match in BCMatch:
        perio = Internal.getNodesFromName(match, 'GridConnectivityProperty')
        if perio != [] :
            rotCenter  = Internal.getNodeFromName(perio, 'RotationCenter')[1]
            rotAngle   = Internal.getNodeFromName(perio, 'RotationAngle')[1]
            translation = Internal.getNodeFromName(perio, 'Translation')[1]

            list1 = [abs(translation).tolist(),abs(rotAngle).tolist(),rotCenter.tolist()]

            if list1 not in liste_BCPeriodiques:

                 liste_BCPeriodiques.append(list1)
                 print(liste_BCPeriodiques)


    t = Internal.rmGhostCells(t, t, 2)
    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')

    for b in bases:
        zones += Internal.getNodesFromType1(b, 'Zone_t')

    somme = 0
    numZone.append(somme)

    for z in zones:
        dim = Internal.getZoneDim(z)
        ni = dim[1]-1
        nj = dim[2]-1
        nk = dim[3]-1

        dimPb = dim[4]

        if dimPb==2:nk=1

        nbbloc_i = ni/taille_bloc
        nbbloc_i = int(nbbloc_i)
        if nbbloc_i*taille_bloc < ni: nbbloc_i = nbbloc_i + 1
        nbbloc_j = nj/taille_bloc
        nbbloc_j = int(nbbloc_j)
        if nbbloc_j*taille_bloc < nj: nbbloc_j = nbbloc_j + 1
        nbbloc_k = nk/taille_bloc
        nbbloc_k = int(nbbloc_k)
        if nbbloc_k*taille_bloc < nk: nbbloc_k = nbbloc_k + 1

        somme += nbbloc_i*nbbloc_j*nbbloc_k

        numZone.append(somme)

        print(nbbloc_i,nbbloc_j,nbbloc_k)

        print('decoupage de la zone', z[0])

        for k in range(0,nbbloc_k):
            for j in range(0,nbbloc_j):
                for i in range(0,nbbloc_i):

                    imin = i*taille_bloc + 1
                    jmin = j*taille_bloc + 1
                    kmin = k*taille_bloc + 1

                    imax = min((i+1)*taille_bloc+1,ni+1)
                    jmax = min((j+1)*taille_bloc+1,nj+1)
                    kmax = min((k+1)*taille_bloc+1,nk+1)

                    if dimPb == 2: kmax=1

                    zp = T.subzone(z,(imin,jmin,kmin),(imax,jmax,kmax))

                    cflmax_loc = C.getMaxValue(zp, 'centers:CFL')
                    dtmin_loc = 1./cflmax_loc
                    if dtmin_loc < dtmin: dtmin = dtmin_loc

                    zones_decoupe.append(zp)

    zones_decoupe_=[]
    for zp in zones_decoupe:
        cflmax_loc = C.getMaxValue(zp, 'centers:CFL')
        dtmin_loc = 1./cflmax_loc
        niveauTps = math.log(((2**exposantMax)*dtmin)/dtmin_loc)/math.log(2)
        if niveauTps < 0.0: niveauTps = 0
        #print(niveauTps, math.ceil(niveauTps))
        niveauTps = math.ceil(niveauTps)
        dicoTps[zp[0]]=niveauTps
        #print('niveau_tps= ', float(pow(2,niveauTps)))
        zp = C._initVars(zp,'centers:niveaux_temps',float(niveauTps))


    print('Nb de zones= ', len(zones_decoupe))


    t[2][1][2] = zones_decoupe

    base=Internal.getBases(t)
    base[0][2].append(fes)
    base[0][2].append(rs)

    t = X.connectMatch(t, tol=1.e-7, dim=dimPb)

    for perio in liste_BCPeriodiques:
        #print('perio[0]= ', perio[0])
        t = X.connectMatchPeriodic(t, rotationCenter=perio[2],rotationAngle=[perio[1][0],perio[1][1],perio[1][2]], translation=perio[0],tol = 1.e-7,dim=dimPb)



    newZones = Internal.getZones(t)


    ### On controle si les zones en contact ont, au plus, un rapport 2 entre leur niveau en temps ###
    ### Tant que ce n'est pas le cas, on modifie les niveaux afin d arriver dans cette situtation ###

    chgt = 1

    while chgt != 0:

        chgt = 0
        for z in newZones:

            gcs = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for gc in gcs:

                name = Internal.getValue(gc)

                if dicoTps[z[0]] - dicoTps[name] > 1:
                    dicoTps[name] = dicoTps[z[0]] - 1
                    node = Internal.getNodeFromName(t, name)
                    C._initVars(node, 'centers:niveaux_temps', dicoTps[name])
                    chgt = 1

            gcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for gc in gcs:

                name = Internal.getValue(gc)

                if dicoTps[z[0]] - dicoTps[name] > 1:
                    dicoTps[name] = dicoTps[z[0]] - 1
                    node = Internal.getNodeFromName(t, name)
                    C._initVars(node, 'centers:niveaux_temps', dicoTps[name])
                    chgt = 1

        print('chgt= ', chgt)

    #### Il se peut que des niveaux les plus bas aient disparu avec les chgts faits pour avor un rapport 2 entre niveaux adjacents #####
    #### On met donc a jour les niveaux en temps ###

    niveau_min = 100000
    for z in newZones:
         niveau = C.getMaxValue(z, 'centers:niveaux_temps')
         if niveau < niveau_min :
             niveau_min = niveau

    print('niveau_min= ', niveau_min)
    if niveau_min != 0:
        for z in newZones:
            niveau_precedent = C.getMaxValue(z, 'centers:niveaux_temps')
            C._initVars(z, 'centers:niveaux_temps', niveau_precedent - niveau_min)



    #### On fusionne lorsque c est possible les zones de meme pas de temps ####

    zones_merged=[]
    for i in range(len(zones)):
        print(numZone[i], numZone[i+1])
        zones_a_meger =  zones_decoupe[numZone[i]:numZone[i+1]]

        niveauMax = 0
        niveauMin = 1000
        dico_niveau={}

        for za in zones_a_meger:

            niveau = C.getMaxValue(za, 'centers:niveaux_temps')

            if niveau > niveauMax : niveauMax = niveau
            if niveau < niveauMin : niveauMin = niveau

        niveauMin = int(niveauMin)
        niveauMax = int(niveauMax)


        for niv in range(niveauMin,niveauMax+1):
            iso_niveau=[]
            for za in zones_a_meger:
                if C.getMaxValue(za, 'centers:niveaux_temps') == niv:
                    iso_niveau.append(za)

            print('nb zones av merge', len(iso_niveau))
            iso_niveau_ = T.merge(iso_niveau, tol=1.e-11)
            print('nb zones ap merge', len(iso_niveau_))

            for ziso in iso_niveau_:
                zones_merged.append(ziso)


    print('nb zones final', len(zones_merged))
    t[2][1][2] = zones_merged


    zones = Internal.getZones(t)
    for z in zones:
        niveau = C.getMaxValue(z, 'centers:niveaux_temps')
        C._initVars(z, 'centers:niveaux_temps', pow(2,niveau))


    if NP > 0 :
        t = _distribMpiDtloc(t, pow(2,niveauMax), NP)


    base=Internal.getBases(t)
    base[0][2].append(fes)
    base[0][2].append(rs)

    t = X.connectMatch(t, tol=1.e-6, dim=dimPb)
    for perio in liste_BCPeriodiques:
        t = X.connectMatchPeriodic(t, rotationCenter=perio[2],rotationAngle=[perio[1][0],perio[1][1],perio[1][2]], translation=perio[0],tol = 1.e-7,dim=dimPb)

    C.convertPyTree2File(t, 'essai.cgns')

    return t

## NOTE: _decoupe4 vs _decoupe2
## _decoupe4 is heavy based on _decoupe2 and is just a "cleaning" of what is done in the latter.
## _decoupe2 is left for now and it will be decided in the future which to keep.
## _decoupe4 uses Cassiopee functions to split the domain and uses the block size algorithm of _decoupe2 to determine
## the number of "cuts" in each direction.
def splitting_per_direction(t,dir,taille_bloc):
    import Transform.PyTree as T
    zones  = Internal.getNodesFromType(t, 'Zone_t')
    zones_delete = []
    for z in zones:
        dim = Internal.getZoneDim(z)
        ni = dim[1]-1
        nj = dim[2]-1
        nk = dim[3]-1

        dimPb = dim[4]
        
        if dimPb == 2: nk=1

        nbbloc_i = ni/taille_bloc
        nbbloc_i = int(nbbloc_i)
        if nbbloc_i*taille_bloc < ni: nbbloc_i = nbbloc_i + 1
        nbbloc_j = nj/taille_bloc
        nbbloc_j = int(nbbloc_j)
        if nbbloc_j*taille_bloc < nj: nbbloc_j = nbbloc_j + 1
        nbbloc_k = nk/taille_bloc
        nbbloc_k = int(nbbloc_k)
        if nbbloc_k*taille_bloc < nk: nbbloc_k = nbbloc_k + 1

        if dir == 1:
            nbloc=nbbloc_i
            nbbloc_j = 1
            nbbloc_k = 1
        elif dir ==2:
            nbbloc_i = 1
            nbloc=nbbloc_j
            nbbloc_k = 1
        else:
            nbbloc_i = 1
            nbbloc_j = 1
            nbloc=nbbloc_k
        if nbloc==1:
            continue
        zones_delete.append(z)
        print('Splitting of Zone %s in num. blocs in i,j,k::%d %d %d'%(z[0],nbbloc_i,nbbloc_j,nbbloc_k))       
        zones_decoupe  = T.splitNParts(z,nbloc, multigrid=0, dirs=[dir],recoverBC=True)
        t[2][1][2] += zones_decoupe
        
    for z in zones_delete:
        Internal._rmNode(t,z)

    return t


def _decoupe4(t,tc=None,exposantMax=2,NP=0,taille_bloc=25,isOctree=False):
    import Transform.PyTree as T

    dimPb = 3
    node = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dimPb = Internal.getValue(node) 

    if not isOctree:
        print("Saving periodic settings...start")
        liste_BCPeriodiques = []
        BCMatch = Internal.getNodesFromType(t, 'GridConnectivity1to1_t')
        for match in BCMatch:
            perio = Internal.getNodesFromName(match, 'GridConnectivityProperty')
            if perio != [] :
                rotCenter  = Internal.getNodeFromName(perio, 'RotationCenter')[1]
                rotAngle   = Internal.getNodeFromName(perio, 'RotationAngle')[1]
                translation = Internal.getNodeFromName(perio, 'Translation')[1]
             
                list1 = [abs(translation).tolist(),abs(rotAngle).tolist(),rotCenter.tolist()]
        
                if list1 not in liste_BCPeriodiques:            
                     liste_BCPeriodiques.append(list1)
                     print("liste_BCPeriodiques=",liste_BCPeriodiques)
        print("Saving periodic settings...end")
    
        t = Internal.rmGhostCells(t, t, 2,adaptBCs=1)
        zones  = Internal.getNodesFromType(t     , 'Zone_t')
        for z in zones:
            gcs = Internal.getNodeFromName(z, 'ZoneGridConnectivity')
            Internal._rmNode(z,gcs)
           
        print("Splitting mesh...start")   
        ## Doing X direction splitting
        dir = 1
        t   = splitting_per_direction(t,dir,taille_bloc)
        ## Doing Y direction splitting
        dir = 2
        t   = splitting_per_direction(t,dir,taille_bloc)
        ## Doing Z direction splitting
        dir = 3
        t   = splitting_per_direction(t,dir,taille_bloc)
        print("Splitting mesh...done")
        
    zones  = Internal.getNodesFromType(t     , 'Zone_t')
    print('Post split: Total number of zones=', len(zones))

    print("Setting 'centers:niveaux_temps' as a flow variables...start")
    dtmin = 1e15
    count = 0
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dimPb == 3:
            zp = T.subzone(z,(3,3,3),(dim[1]-2,dim[2]-2,dim[3]-2))
        else:
            zp = T.subzone(z,(3,3,1),(dim[1]-2,dim[2]-2,1))
        cflmax_loc = C.getMaxValue(zp, 'centers:CFL')
        dtmin_loc = 1./cflmax_loc
        if dtmin_loc < dtmin: dtmin = dtmin_loc

    if not isOctree:
        for z in zones:
            gcs = Internal.getNodeFromName(z, 'ZoneGridConnectivity')
            Internal._rmNode(z,gcs)
        
        t = X.connectMatch(t, tol=1.e-7, dim=dimPb)
        for perio in liste_BCPeriodiques:
            t = X.connectMatchPeriodic(t, rotationCenter=perio[2],rotationAngle=[perio[1][0],perio[1][1],perio[1][2]], translation=perio[0],tol = 1.e-7,dim=dimPb)
   
    dicoTps = {}
    zones  = Internal.getNodesFromType(t, 'Zone_t')
    ## Note: 2 parameters to control the number of time zones
    ## 1) taille_bloc: The smaller this number is the more likely you will have the same
    ##                  number of time levels as exposantMax
    ## 2) exposantMax: for forcing niveauTps < 0.0 and by making exposantMax< thexposantMax
    ##                 one can play with taille_bloc such that we have more control of the regions of
    ##                 with smaller time steps without "exploding" the different levels of time zones.
    for z in zones:
        cflmax_loc = C.getMaxValue(z, 'centers:CFL')
        dtmin_loc = 1./cflmax_loc
        niveauTps = math.log2(((2**exposantMax)*dtmin)/dtmin_loc)
        if niveauTps < 0.0: niveauTps = 0
        niveauTps = math.ceil(niveauTps)
        dicoTps[z[0]]=niveauTps
        z = C._initVars(z,'centers:niveaux_temps',float(niveauTps))
        
    for z in zones:
        cflmax_loc = C.getMaxValue(z, 'centers:CFL')
        dtmin_loc = 1./cflmax_loc
        niveauTps = math.log2(((2**exposantMax)*dtmin)/dtmin_loc)    
        z = C._initVars(z,'centers:niveaux_tempsP1',niveauTps)
        
    print("Setting 'centers:niveaux_temps' as a flow variables...done")

    ### On controle si les zones en contact ont, au plus, un rapport 2 entre leur niveau en temps ###
    ### Tant que ce n'est pas le cas, on modifie les niveaux afin d arriver dans cette situtation ###
    print("Checking smooth transition in times levels...start")
    chgt = 1    
    while chgt != 0:
        chgt = 0      
        for z in zones:
            gcs = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for gc in gcs:
                name = Internal.getValue(gc)                
                if dicoTps[z[0]] - dicoTps[name] > 1:
                    dicoTps[name] = dicoTps[z[0]] - 1
                    node = Internal.getNodeFromName(t, name)
                    C._initVars(node, 'centers:niveaux_temps', dicoTps[name])
                    chgt = 1

            gcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for gc in gcs: 
                name = Internal.getValue(gc)
                if dicoTps[z[0]] - dicoTps[name] > 1:
                    dicoTps[name] = dicoTps[z[0]] - 1
                    node = Internal.getNodeFromName(t, name)
                    C._initVars(node, 'centers:niveaux_temps', dicoTps[name])
                    chgt = 1
           
    print("Checking smooth transition in times levels...done")
    
    niveau_min = 100000
    max_time_level=0
    for z in zones:
         niveau = C.getMaxValue(z, 'centers:niveaux_temps')
         if niveau < niveau_min:
             niveau_min = niveau
    print('Smallest time level= ', niveau_min)
    
    print("Setting smallest time level to 0...start")
    if niveau_min != 0:
        for z in zones:
            niveau_precedent = C.getMaxValue(z, 'centers:niveaux_temps')         
            C._initVars(z, 'centers:niveaux_temps', niveau_precedent - niveau_min)

    for z in zones:
        niveau_precedent = C.getMaxValue(z, 'centers:niveaux_temps')
        if niveau_precedent>max_time_level: max_time_level=niveau_precedent
    print("Setting smallest time level to 0...done")
    
    print("Total number of time levels=",max_time_level)
    if not isOctree:
        print("Number of zones pre merge=",len(zones))
        
        print("Merging zones with same time level...start")
        count = 0
        for i in range(0,int(max_time_level)+1):
            list_same_levels=[]
            for z in zones:
                current_level = C.getMaxValue(z, 'centers:niveaux_temps')   
                if current_level == i:
                    list_same_levels.append(z)
            new_zone=T.merge(list_same_levels)
            t[2][1][2] += new_zone
            count += 1
            del list_same_levels
            
        for z in zones:
            Internal._rmNode(t,z)
            
        zones  = Internal.getNodesFromType(t, 'Zone_t')
        print("Merging zones with same time level...done")
        print("Number of zones post merge=",len(zones))
        print("Dimension=",dimPb)

    for z in zones:
        niveau = C.getMaxValue(z, 'centers:niveaux_temps')
        C._initVars(z, 'centers:niveaux_temps', pow(2,niveau))


    if not isOctree:
        ##This has not been tested
        if NP > 0 :
            t = _distribMpiDtloc(t, pow(2,niveauMax), NP)
    
        zones = Internal.getNodesFromType(t, 'Zone_t')
        for z in zones:
            gcs = Internal.getNodeFromName(z, 'ZoneGridConnectivity')
            Internal._rmNode(z,gcs)
            
        t = X.connectMatch(t, tol=1.e-7, dim=dimPb)
        for perio in liste_BCPeriodiques:
            t = X.connectMatchPeriodic(t, rotationCenter=perio[2],rotationAngle=[perio[1][0],perio[1][1],perio[1][2]], translation=perio[0],tol = 1.e-7, dim=dimPb)
    
        C.addState2Node__(t, 'EquationDimension', dimPb)   
        t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
        if dim == 2: 
            t = T.addkplane(t)
            t = T.contract(t, (0,0,0), (1,0,0), (0,1,0),0.025)
            t = T.makeDirect(t)
        tc = C.node2Center(t)
        tc = X.setInterpData2(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, method='lagrangian',dim=dim)
        tc = C.rmVars(tc, 'FlowSolution')
        tc = C.rmVars(tc, 'CellN')
    else:
        tc_orig = Internal.copyRef(tc)
        tc = C.node2Center(t)
        tc = X.setInterpData2(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, method='lagrangian',dim=dim)
        bases = Internal.getNodesFromType(tc_orig, 'CGNSBase_t')
        for b in bases:
            zones = Internal.getNodesFromType(b, 'Zone_t')
            for z in zones:
                    zones2 = Internal.getNodesFromType(z, 'ZoneSubRegion_t')
                    for z2 in zones2:
                        if z2[0][0:3]=='IBC':
                            node = Internal.createNode('PointRange', 'IndexArray_t', value=[0,0,0,0,0,0], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('PointRangeDonor', 'IndexArray_t', value=[0,0,0,0,0,0], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('DirDonneur', 'IndexArray_t', value=[0], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('DirReceveur', 'IndexArray_t', value=[0], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('Transform', 'IndexArray_t', value=[1,2,3], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('PointPivot', 'IndexArray_t', value=[0,0,0], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('Profondeur', 'IndexArray_t', value=[2], children=[])
                            Internal.addChild(z2, node, pos=-1)                
                            node = Internal.createNode('NMratio', 'IndexArray_t', value=[1,1,1], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('DnrZoneName', 'IndexArray_t', value=[Internal.getValue(z2)], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            
                            t_zone     = Internal.getNodesFromName(t, Internal.getValue(z2))
                            time_level = C.getMaxValue(t_zone, 'centers:niveaux_temps')
            
                            node = Internal.createNode('LevelZRcv', 'IndexArray_t', value=[time_level], children=[])
                            Internal.addChild(z2, node, pos=-1)
                            node = Internal.createNode('LevelZDnr', 'IndexArray_t', value=[time_level], children=[])
                            Internal.addChild(z2, node, pos=-1)
            
                            path = b[0]+'/'+z[0]
                            node = Internal.getNodeFromPath(tc, path)
                            Internal.addChild(node, z2, pos=-1)        
    return (t,tc)


def set_dt_lts(t,running_cfl=None,dt=None):
    import Transform.PyTree as T
    if running_cfl is None and dt is None:
        raise ValueError("set_dts_lts: running_cfl and dt cannot be both None. Either running_cfl or dt for the most refined level must be provided.")
        
    dimPb = 3
    node  = Internal.getNodeFromName(t, 'EquationDimension')
    if node is not None: dimPb = Internal.getValue(node) 
    time_level_glob = 0
    cflmax_glob     = 0.0
    zones           = Internal.getZones(t)
    for z in zones:
        dim        = Internal.getZoneDim(z)
        if dimPb == 3:
            z          = T.subzone(z, (3,3,3), (dim[1]-2,dim[2]-2,dim[3]-2))
        else:
            z          = T.subzone(z, (3,3,1), (dim[1]-2,dim[2]-2,1))
            
        cflmax     = C.getMaxValue(z, 'centers:CFL')
        time_level = C.getMaxValue(z, 'centers:niveaux_temps')
        if time_level > time_level_glob:time_level_glob = time_level
        if cflmax > cflmax_glob:cflmax_glob = cflmax
    
    dt_max = running_cfl*(1/cflmax_glob)
    if dt is not None:
        dt_max = dt
    print("dt_refined=",dt_max)
    dt_max = dt_max*2**(numpy.log2(time_level_glob))
    print("Total number of time levels=",time_level_glob)
    print("dt_coarse=",dt_max)
    
    return dt_max

#==============================================================================
# distribution mpi pour dtloc instationnaire
#==============================================================================
def _distribMpiDtloc(t,niveauMax,NP):

    import Transform.PyTree as T
    import Distributor2.PyTree as D2

    list_level= []

    zones = Internal.getZones(t)

    exposant_max = math.log(niveauMax)/math.log(2)
    exposant_max = int(exposant_max)

    for exposant in range(0,exposant_max+1):
      list=[]
      for z in zones:
        dim = Internal.getZoneDim(z)
        i = dim[1]//2
        j = dim[2]//2
        k = dim[3]//2
        level = C.getValue( z, 'centers:niveaux_temps', (i,j,k))
        if pow(2,exposant) == int(level):
             list.append(z) ### list contient toutes les zones de meme niveau en temps
      list_level.append(list)

    tp = C.newPyTree(['Base'])

    for level in list_level: ### On distribue chaque niveau en tps de maniere independante sur NP proc
      level = T.splitSize(level,N=0, multigrid=0, dirs=[1,2,3], R=NP, minPtsPerDir=5)
      level = X.connectMatch(level)
      level, stats = D2.distribute(level, NP,algorithm='gradient')

      print(stats)

      tp[2][1][2] += level


    return tp

#==============================================================================
# compute_dpJ_dpW in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute_dpJ_dpW(t, teff, metrics, cosAoA, sinAoA, surfinv):

    zones     = Internal.getZones(t)
    zones_eff = Internal.getZones(teff)

    if FastC.HOOK is None:
            own    = Internal.getNodeFromName1(t   , '.Solver#ownData')  # noeud
            dtloc  = Internal.getNodeFromName1(own , '.Solver#dtloc')    # noeud
            dtloc  = Internal.getValue(dtloc)                       # tab numpy
            FastC.HOOK   = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)
            nitrun =0; nstep =1
            hook1  = FastC.HOOK.copy()
            hook1.update(FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, 0) )
    else:   hook1  = FastC.HOOK

    fasts.compute_dpJ_dpW(zones, zones_eff, metrics, hook1, cosAoA, sinAoA, surfinv)
    return None


'''
#==============================================================================
# computeAdjoint in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _computeAdjoint(t, metrics, nit_adjoint, indFunc, tc=None, graph=None):

    bases  = Internal.getNodesFromType1(t  , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(t    , '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own  , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t')

    node = Internal.getNodeFromName1(   t, '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = FastC.OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    dtloc = Internal.getValue(dtloc) # tab numpy
    nitmax = int(dtloc[0])


    if tc is not None:
         bases = Internal.getNodesFromType1(tc, 'CGNSBase_t')  # noeud
         tc_compact = Internal.getNodeFromName1(bases[0], 'Parameter_real')
         if tc_compact is not None:
                param_real= tc_compact[1]
                param_int = Internal.getNodeFromName1( bases[0], 'Parameter_int' )[1]

                zonesD = []
                for f in bases:
                    tmp = Internal.getNodesFromType1(f, 'Zone_t')
                    zonesD += tmp


    return None

#==============================================================================
# computedJdX in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _computedJdX(t, metrics, nitrun, tc=None, graph=None, indFunc):


    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t')

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    omp_mode = 0
    if  node is not None: omp_mode = Internal.getValue(node)

    dtloc = Internal.getValue(dtloc) # tab numpy
    nitmax = int(dtloc[0])


    if tc is not None:
         bases = Internal.getNodesFromType1(tc, 'CGNSBase_t')  # noeud
         tc_compact = Internal.getNodeFromName1(bases[0], 'Parameter_real')
         if tc_compact is not None:
                param_real= tc_compact[1]
                param_int = Internal.getNodeFromName1( bases[0], 'Parameter_int' )[1]

                zonesD = []
                for f in bases:
                    tmp = Internal.getNodesFromType1(f, 'Zone_t')
                    zonesD += tmp

   #--------------------------------------------------------------------------------

    fasts._computedJdX(zones, metrics, omp_mode, hook1, indFunc)

    return None

'''


##ADD to the t_conv_history
def add2previous(var_list,zone_save,zone_current,it0):
    
    new_value_node = Internal.getValue(zone_current) + Internal.getValue(zone_save)
    
    for var_local in var_list:
        RSD_tmp   =Internal.getNodeByName(zone_current,var_local)[1]
        RSD_sv    =Internal.getNodeByName(zone_save,var_local)[1]

        sh_local  = numpy.shape(RSD_sv)[0]

        nrec        = Internal.getValue(zone_save)
        nrec_tmp    = Internal.getValue(zone_current)
               
        sh_tmp = nrec_tmp*sh_local//nrec
        sh     = nrec    *sh_local//nrec        
        
        local_val = numpy.zeros(sh+sh_tmp)
    
        shift_local=0
        if var_local == 'IterationNumber': shift_local = it0
        for i in range(sh): local_val[i]=RSD_sv[i]
        for i in range(sh_tmp): local_val[i+sh]=RSD_tmp[i]+shift_local    
        Internal.createUniqueChild(zone_save,var_local,'DataArray_t',value=local_val)
    Internal.setValue(zone_save,new_value_node)
    return None


##POPULATE GLOBAL CONVERGENCE FOR EACH BASE
def calc_global_convergence(t):
    var_list    =['RSD_oo','RSD_L2','RSD_oo_diff','RSD_L2_diff']
    var_list_L2 =['RSD_L2','RSD_L2_diff']
    var_list_Loo=['RSD_oo','RSD_oo_diff']
    neq=5
    if Internal.getValue(Internal.getNodeFromType(t, 'GoverningEquations_t'))=='NSTurbulent':
        neq = 6

    for b in Internal.getBases(t):
        nzones= len(Internal.getZones(b))
        c     = Internal.getNodeByName(b,'GlobalConvergenceHistory')
        
        
        zones                 = Internal.getZones(b)
        zone_cong_history     = Internal.getNodeByName(zones[0],'ZoneConvergenceHistory')
        RSD_                  = Internal.getNodeByName(zone_cong_history,'IterationNumber')[1]
        nrec                  = Internal.getValue(zone_cong_history)

        ##NOTE: the value of the ZoneConvergenceHistory node has the number of rec taking at time t
        ##      the shape of the array of IterationNumber is total number of nrec.
                
        Internal.setValue(c, nrec);

        ##IterationNumber Node
        Internal.createChild(c, 'IterationNumber', 'DataArray_t', numpy.zeros((nrec), numpy.int32))
        RSD_b                 = Internal.getNodeByName(c,'IterationNumber')[1]
        for i in range(nrec):RSD_b[i]=RSD_[i]

        ##The rest of the nodes in the convergence history
        for var in var_list:
            tmp = numpy.zeros((nrec*neq), numpy.float64)
            Internal.createChild(c, var ,'DataArray_t', tmp)                    

        total_Ncells = 0
        for z in Internal.getZones(b):
            #if z[0]=='GlobalConvergentHistory': continue
            zone_cong_history     =Internal.getNodeByName(z,'ZoneConvergenceHistory')            
            isCalcCheck           =Internal.getNodeByName(z,'isCalc')
            
            if isCalcCheck is None:isCalc=1
            else: isCalc=Internal.getValue(isCalcCheck)

            nx = Internal.getValue(z)[0][0]
            ny = Internal.getValue(z)[1][0]
            nz = Internal.getValue(z)[2][0]
            total_Ncells += nx*ny*nz*isCalc

            if isCalc==1:
                ##Loo per base
                for var_local in var_list_Loo:
                    RSD_     = Internal.getNodeByName(zone_cong_history,var_local)[1]
                    RSD_b    = Internal.getNodeByName(c,var_local)[1]
                    sh_local = numpy.shape(RSD_b)[0]
                    sh       = nrec    *sh_local//nrec
                    for i in range(sh):
                        RSD_b[i]=max(RSD_b[i],RSD_[i])
                
                ##L2 per base (here only sum)
                for var_local in var_list_L2:
                    RSD_     = Internal.getNodeByName(zone_cong_history,var_local)[1]
                    RSD_b    = Internal.getNodeByName(c,var_local)[1]
                    sh_local = numpy.shape(RSD_b)[0]
                    sh       = nrec    *sh_local//nrec
                    for i in range(sh):
                        RSD_b[i]+=RSD_[i]*nx*ny*nz

        ##Average L2 per base            
        for var_local in var_list_L2:
                RSD_b    = Internal.getNodeByName(c,var_local)[1]
                RSD_b[:] = RSD_b[:]/total_Ncells #nzones            
    return t


##ADD or CREATE a t_conv_history (tree that stores only the convergence history)
def create_add_t_converg_hist(t2,it0=0,t_conv_hist=None):
    var_list    =['IterationNumber','RSD_oo','RSD_L2','RSD_oo_diff','RSD_L2_diff']
    var_list_L2 =['RSD_L2','RSD_L2_diff']
    var_list_Loo=['RSD_oo','RSD_oo_diff']
    
    ##Create local copy by removing unnecessary stuff
    t=Internal.rmNodesByName(t2,'ReferenceState')
    list_vars=['.Solver#define','.Solver#ownData']
    for var in list_vars:Internal._rmNodesByName(t,var)

    ##Remove unncessary nodes in the zones    
    for z in Internal.getZones(t):
        listzones_rm=[]

        CalcNode = Internal.getNodeFromName(z, 'isCalc')
        if CalcNode is None:
            ##CHECK to see if cellN==1 can be found in the zone
            ## if N_cells(cellN=1)=0: isCalc=0
            ## else: isCalc=1
            sol = Internal.getNodeFromName(z, 'FlowSolution#Centers')
            if sol:
                celN           = Internal.getNodeFromName(sol, 'cellN')[1]
                cells_zone_1   = numpy.count_nonzero(celN==1)
                isCalc = 0
                if cells_zone_1>0: isCalc=1
            else:
                isCalc = 1
            Internal.createUniqueChild(z, 'isCalc', 'DataArray_t', value=isCalc)
        
        for z2 in Internal.getChildren(z):
            if z2[0] != 'ZoneConvergenceHistory' and z2[0] != 'ZoneType' and z2[0] != 'isCalc':
                listzones_rm.append(z2)
        for z2 in listzones_rm:
            Internal._rmNode(z,z2)
       
    ##Check to see if I need to write the file or append    
    if t_conv_hist is None:
        if it0>0:
            shift_local    = it0
            for z in Internal.getZones(t):
                RSD       = Internal.getNodeByName(z,'IterationNumber')[1]
                sh_local  = numpy.shape(RSD)[0]                
                for i in range(sh_local):
                    RSD[i]=RSD[i]+shift_local
            RSD       = Internal.getNodeByName(Internal.getNodeByName(t,'GlobalConvergenceHistory'),'IterationNumber')[1]
            sh_local  = numpy.shape(RSD)[0]                
            for i in range(sh_local):
                RSD[i]=RSD[i]+shift_local
        return t    
    else:
        for z in Internal.getZones(t_conv_hist):
            zone_save     =Internal.getNodeByName(z,'ZoneConvergenceHistory')
            zone_current  =Internal.getNodeByName(Internal.getNodeByName(t,z[0]),'ZoneConvergenceHistory')
            add2previous(var_list,zone_save,zone_current,it0)
            
        for z in Internal.getBases(t_conv_hist):
            zone_save    = Internal.getNodeByName(z,'GlobalConvergenceHistory')
            
            zone_current = Internal.getNodeByName(Internal.getNodeByName(t,z[0]),'GlobalConvergenceHistory')
            if Internal.getChildren(zone_current) and Internal.getChildren(zone_save):
                add2previous(var_list,zone_save,zone_current,it0)
            
        return t_conv_hist


