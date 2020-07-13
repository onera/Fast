"""Fast Structured Grid Navier-Stokes solver.
"""
from . import fasts
from . import FastS
__version__ = FastS.__version__

OMP_MODE = 0

import numpy
import os
try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Connector
    import Connector.PyTree as X
    import FastC.PyTree as FastC
    import KCore
    import math
    import time as Time

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
def _stretch(coord,nbp,x1,x2,dx1,dx2,ityp):
    fasts._stretch(coord,nbp,x1,x2,dx1,dx2,ityp)
    return None

#==============================================================================
# compute in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c", NIT=1):
#def _compute(t, metrics, nitrun, tc=None, graph=None, layer="Python", NIT=1):
    """Compute a given number of iterations."""

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t') 

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    omp_node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode  = OMP_MODE
    if  omp_node is not None: ompmode = Internal.getValue(omp_node)

    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])                 

    #### a blinder...
    param_int_firstZone = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1]
    itypcp = param_int_firstZone[29]
    rk     = param_int_firstZone[52]
    exploc = param_int_firstZone[54]
    #### a blinder...
                
     
    if nitrun == 1: print('Info: using layer trans=%s (ompmode=%d)'%(layer, ompmode))

    if layer == "Python": 
        
      if exploc==2 and tc is not None and rk==3:
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
       
      for nstep in range(1, nitmax+1): # pas RK ou ssiterations

         #print('nstep=%d'%nstep)

         hook1 = FastC.HOOK.copy()
         distrib_omp = 0
         hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )
         nidom_loc = hook1["nidom_tot"]
        
         skip      = 0
         if hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1: skip = 1

         # calcul Navier_stokes + appli CL
         if nidom_loc > 0 and skip == 0:
            nstep_deb = nstep
            nstep_fin = nstep
            layer_mode= 0
            nit_c     = 1
            #t0=Time.time()
            fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, hook1)
            #print('t_compute = %f'%(Time.time() - t0))

            # dtloc GJeanmasson
            if exploc==2 and tc is not None and rk==3:
               fasts.dtlocal2para_(zones, zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1, dest)

               if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
               elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
               _applyBC(zones,metrics, hook1, nstep, ompmode, var=vars[0])    

               FastC.switchPointers2__(zones,nitmax,nstep)
                
               # Ghostcell
               if nitmax%3 != 0: # Tous les schemas sauf constantinescu RK3
                   if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
                   elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1']
                   timelevel_target = int(dtloc[4])
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, nitmax, rk, exploc)
            
               fasts.recup3para_(zones,zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1) 

               if nstep%2 == 0:
                   timelevel_target = int(dtloc[4])
                   vars = ['Density'  ]
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, nitmax, hook1, nitmax, rk, exploc, 2)

               if    nstep%2 == 0 and itypcp == 2 : vars = ['Density'  ] 
               elif  nstep%2 == 1 and itypcp == 2 : vars = ['Density_P1'] 
               _applyBC(zones,metrics, hook1, nstep, ompmode,  var=vars[0])



            else:
              #Ghostcell
              vars = FastC.varsP
              if nstep%2 == 0 and itypcp == 2 : vars = FastC.varsN  # Choix du tableau pour application transfer et BC
              timelevel_target = int(dtloc[4])
              #t0=Time.time()
              _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1)
              #print('t_transferts = %f'%(Time.time() - t0)

      dtloc[3] +=1    #time_level_motion
      dtloc[4] +=1    #time_level_target
    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT

      if tc is None:
          FastC.HOOK['param_int_tc']  = None
          FastC.HOOK['param_real_tc'] = None

      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c , FastC.HOOK)


    #switch pointer a la fin du pas de temps
    if exploc==2 and tc is not None and rk==3:
         if layer == 'Python':
             FastC.switchPointers__(zones, 1, 3)
         else:
             FastC.switchPointers3__(zones,nitmax)
    else:
         case = NIT%3
         if case != 0 and itypcp < 2: FastC.switchPointers__(zones, case)
         if case != 0 and itypcp ==2: FastC.switchPointers__(zones, case, nitmax%2+2)

    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT          = 1

    return None


#==============================================================================
# compute in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute_matvec(t, metrics, no_vect_test, tc=None, graph=None):
    

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t') 

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    nstep = 1

    hook1       = FastC.HOOK.copy()
    distrib_omp = 0
    nitrun      = 0
    hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )
        
    fasts._matvecPT(zones, metrics, nitrun, no_vect_test, ompmode, hook1)

    return None
#
#==============================================================================
# alloue retourne la metrique
#==============================================================================
def allocate_metric(t):
    zones        = Internal.getZones(t)
    dtloc        = Internal.getNodeFromName3(t, '.Solver#dtloc')
    dtloc_numpy  = Internal.getValue(dtloc)
    nssiter      = int(dtloc_numpy[0])

    metrics=[]; motion ='none'
    for z in zones:
        b = Internal.getNodeFromName2(z, 'motion')
        if b is not None: motion = Internal.getValue(b)
        num = Internal.getNodeFromName1(z, '.Solver#ownData')
        if num is None:
            raise ValueError("metric: numerics is missing for zone %s."%z[0])
        if motion == 'rigid' or  motion == 'deformation':
            grids = Internal.getNodesFromType1(z, 'GridCoordinates_t')
            if len(grids) == 1:
               grid_init = Internal.copyTree(grids[0])
               grid_init[0] = 'GridCoordinates#Init'
               Internal.addChild(z, grid_init, pos=-1) # last
        metrics.append(fasts.allocate_metric(z, nssiter))
    return metrics
#
#==============================================================================
# alloue retourne tableau ssor
#==============================================================================
def allocate_ssor(t, metrics, hook, ompmode):
    zones = Internal.getZones(t)
    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')
    dtloc_numpy = Internal.getValue(dtloc)
    nssiter = int(dtloc_numpy[0])

    ssors = []
    ssors = fasts.allocate_ssor(zones, metrics, nssiter, hook, ompmode)

    return ssors

#
#==============================================================================
# alloue retourne la metrique
#==============================================================================
def _init_metric(t, metrics, omp_mode):
    zones        = Internal.getZones(t)

    fasts.init_metric(zones, metrics, omp_mode)
    
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
def warmup(t, tc, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None):
    """Perform necessary operations for the solver to run."""

    
    # Get omp_mode
    ompmode = OMP_MODE
    node = Internal.getNodeFromName2(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)
        

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastC._reorder(t, tc, ompmode)

    # Construction param_int et param_real des zones
    FastC._buildOwnData(t, Padding)

    #init hook necessaire pour info omp
    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy
    zones = Internal.getZones(t)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it

    # allocation d espace dans param_int pour stockage info openmp
    FastC._build_omp(t) 


    # alloue metric: tijk, ventijk, ssiter_loc
    # init         : ssiter_loc
    metrics = allocate_metric(t)

    # Contruction BC_int et BC_real pour CL
    FastC._BCcompact(t) 


    #determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)
    #print 'int(dtloc[0])= ', int(dtloc[0])
    for nstep in range(1, int(dtloc[0])+1):
        hook1       = FastC.HOOK.copy()
        distrib_omp = 1
        hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, 1, nstep, distrib_omp) )   
    
    #init metric
    #print(ompmode)
    _init_metric(t, metrics, ompmode)

    ssors = allocate_ssor(t, metrics, hook1, ompmode)

    # compact + align + init numa
    rmConsVars=True
    adjoint=Adjoint

    t, FastC.FIRST_IT, zones2compact = FastC.createPrimVars(t, ompmode, rmConsVars, adjoint)
    FastC.HOOK['FIRST_IT']= FastC.FIRST_IT
    #compactage des champs en fonction option de calcul  
    count = -1
    if ompmode == 1: count = 0          
    for data in zones2compact:
        if ompmode == 1: count += 1
        zone    = data[0]
        varnames= data[1]
        for fields in varnames:
            _compact(zone, fields=fields, mode=count)

    #corection pointeur ventijk si ale=0: pointeur Ro perdu par compact.
    zones = Internal.getZones(t) # car create primvar rend zones caduc
    c   = 0
    ale = False
    for z in zones:
        motion = 'none'
        b = Internal.getNodeFromName2(z, 'motion')
        if b is not None: motion = Internal.getValue(b)
        if motion == 'none':
            sol = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
            ro = Internal.getNodeFromName1(sol, 'Density')
            metrics[c][2] = ro[1]
        else: ale = True
        c += 1

    #
    # mise a jour vitesse entrainememnt
    #
    #t0=timeit.default_timer()
    if ale == True and infos_ale is not None:
        print("ale actif. Teta et tetap=", infos_ale)
        teta = infos_ale[0];  tetap = infos_ale[1]
        FastC._motionlaw(t, teta, tetap)
        _computeVelocityAle(t,metrics)
    #t1=timeit.default_timer()
    #print "cout mise a jour vitesse entr= ", t1-t0
 
    #
    # Compactage arbre transfert et mise a jour FastC.HOOK
    #
    if tc is not None:
       X.miseAPlatDonorTree__(zones, tc, graph=graph,list_graph=list_graph) 
    
       FastC.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
       param_real_tc        = Internal.getNodeFromName1( tc, 'Parameter_real')
       if param_real_tc is not None: FastC.HOOK['param_real_tc']= param_real_tc[1]


       # Add linelets arrays in the FastC.HOOK for ODE-based WallModel (IBC)
       nbpts_linelets = 45
       _createTBLESA(tc,nbpts_linelets)
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
        _compact(tmy, fields=varmy)


    #
    # remplissage ghostcells
    #
    hook1["lexit_lu"]= 0
    nstep            = 1
    nitrun           = 0

    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4]) 

    _fillGhostcells(zones, tc, metrics, timelevel_target, ['Density'], nstep,  ompmode, hook1) 
    if tc is not None: C._rmVars(tc, 'FlowSolution')

    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    fasts._computePT_mut(zones, metrics, hook1)

    return (t, tc, metrics)

#==============================================================================
def _compact(t, containers=[Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__], fields=None, mode=None, init=True):
    #global  CACHELINE
    if  mode is not None:
      if mode == -1: thread_numa = 0
      else: thread_numa = 1
    else: thread_numa = 0

    zones = Internal.getZones(t)
    for z in zones:
        ars = FastC.getFields2Compact__(z, containers, fields)
        sh = None ; size = None
        val = [] # valid fields
        for a in ars:
            a1 = a[1]
            if sh is None: sh = a1.shape; size = a1.size; val.append(a)
            elif a1.shape == sh: val.append(a)
        nfields = len(val)
        if nfields > 0:
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')  # noeud
            # Create an equivalent contiguous numpy [flat]
    	    #eq = KCore.empty(size*nfields, CACHELINE)
            eq = numpy.empty(nfields*(size+ param_int[1][66]), dtype=numpy.float64) # add a shift  between prim. variables (param_int[1][SHIFTVAR])
            c = 0        
            if param_int is None:
                raise ValueError("_compact: Parameter_int is missing for zone %s."%z[0])
            for a in val:
                #a[1] = a[1].reshape((size), order='F')
                a1 = a[1]
                ## marc a1 = a[1]
                #print 'a0',a[0],a[1].shape
                # Copy elements
                ptr = a1.reshape((size), order='F') # no copy I hope
                if init: fasts.initNuma( ptr, eq, param_int, c, thread_numa )
                #fasts.initNuma( ptr, eq, param_int, c )
                ## marc ptr = a1.reshape((size), order='F') # no copy I hope
                ## marc fasts.initNuma( ptr, eq, param_int, c )
                #fasts.initNuma( a[1], eq, param_int, c )
                #fasts.initNuma( ptr , eq, param_int, c )
                #eq[c*size:(c+1)*size] = ptr[:]   
                # Replace numpys with slice
                a[1] = eq[c*(size)+c*param_int[1][66]:(c+1)*(size)+c*param_int[1][66]]
                a[1] = a[1].reshape(sh, order='F')
                ## marc a[1] = eq[c*size:(c+1)*size]
                ## marc a[1] = a[1].reshape(sh, order='F')
                #print 'a1',a[0],a[1].shape
                ##a[1] = a[1].reshape( (size), order='F')
                #print 'a',a[0],a[1].shape 

                c += 1
    return None

#==============================================================================
# Prepare IBC ODE (TBLE + Spalart 1D)
#==============================================================================
def _createTBLESA(tc,nbpts_linelets=45):
      
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

                          
         shift        = 0
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
                                       
                                         for noind in range(0,Nbpts_D):
                                             
                                             b0 = valxI[noind]-valxW[noind]
                                             b1 = valyI[noind]-valyW[noind]
                                             b2 = valzI[noind]-valzW[noind]
                                             normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)
                                             
                                             _stretch(linelets[shift + noind*nbpts_linelets: shift + (noind+1)*nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)
                                     else:
   
                                             b0 = valxI-valxW
                                             b1 = valyI-valyW
                                             b2 = valzI-valzW
                                             normb = numpy.sqrt(b0*b0+b1*b1+b2*b2)                                  
                                             
                                             # _stretch(linelets[shift : shift + nbpts_linelets],nbpts_linelets,1.0e-6,0.0,normb)  
                                             _stretch(linelets[shift : shift + nbpts_linelets],nbpts_linelets,0.0,normb, 1.0e-6, 0.6*normb, 2)  
   
   
                                     Internal.createUniqueChild(s, 'CoordinateN_ODE', 'DataArray_t', 
                                                                linelets[shift : shift + Nbpts_D*nbpts_linelets])
                                     Internal.createUniqueChild(s, 'VelocityT_ODE', 'DataArray_t', 
                                                                linelets[shift + Nbpts_D*nbpts_linelets : shift + 2*Nbpts_D*nbpts_linelets ]) 
                                     Internal.createUniqueChild(s, 'TurbulentSANuTilde_ODE', 'DataArray_t', 
                                                                linelets[shift + 2*Nbpts_D*nbpts_linelets  : shift + 3*Nbpts_D*nbpts_linelets ]) 
   
                                        
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
# For periodic unsteady chimera join, parameter must be updated peridicaly 
#==============================================================================
def _UpdateUnsteadyJoinParam(t, tc, omega, timelevelInfos, graph, tc_steady='tc_steady.cgns', directory='.'):

    bases = Internal.getNodesFromType1(t      , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = None
    if own is not None: dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    #on cree les noeud infos insta pour chimere perio s'il n'existe pas 
    TimeLevelOpts=['TimeLevelMotion','TimeLevelTarget'] 
    for opt in TimeLevelOpts:
       tmp = Internal.getNodeFromName1(t, opt)
       if tmp is None: Internal.createUniqueChild(t, opt, 'DataArray_t', value=0)

    if dtloc is not None:
      dtloc            = Internal.getValue(dtloc) # tab numpy
      timelevel_motion = dtloc[3]
      timelevel_target = dtloc[4]
    else:
      timelevel_motion = Internal.getNodeFromName1(t, 'TimeLevelMotion')[1][0]
      timelevel_target = Internal.getNodeFromName1(t, 'TimeLevelTarget')[1][0]

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

       FILE = tc_steady
       if os.access(FILE, os.F_OK): tc  = C.convertFile2PyTree(FILE)
       else: print("error reading %s."%FILE)
       FILE = directory+'/tc_'+str(root)+'.cgns'
       print('File inst=', FILE, 'target=', timelevel_target, 'motionlevel=', timelevel_motion)
       if os.access(FILE, os.F_OK): tc_inst  = C.convertFile2PyTree(FILE)
       else: print("error reading %s."%FILE)

       #
       #timelevel_motion larger than calculated peridicity; need to modify angle of rotation for azymuth periodicity
       #
       if timelevel_motion >= timelevel_period: 
          bases  = Internal.getNodesFromType1(tc_inst , 'CGNSBase_t')       # noeud

          sign =-1
          if omega > 0: sign = 1
          for base in bases:
            if   base[0]=='Rotor': teta = -2*math.pi*timelevel_period/timelevel_360*No_period*sign
            elif base[0]=='Stator':teta =  2*math.pi*timelevel_period/timelevel_360*No_period*sign
            zones  = Internal.getNodesFromType1(base , 'Zone_t')       # noeud
            for z in zones:
              angles = Internal.getNodesFromName2(z, 'RotationAngle')
              for angle in angles: angle[1][:]= angle[1][:] + teta*timelevel_axeRot[:]

       tc = Internal.merge( [tc, tc_inst] )

       # Get omp_mode
       ompmode = OMP_MODE
       node = Internal.getNodeFromName2(t, '.Solver#define')
       if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

       # Reordone les zones pour garantir meme ordre entre t et tc
       FastC._reorder(t, tc, ompmode)

       # Compactage arbre transfert
       g = None; l = None
       zones=Internal.getZones(t)
       X.miseAPlatDonorTree__(zones, tc, procDict=g, procList=l)

       #Remise zero target
       if dtloc is not None: dtloc[4] = 0

    #
    #timelevel_motion larger than number of timelevels for 360degre 
    #
    if timelevel_motion > timelevel_360: dtloc[3] = 0  # remise a zero du compteur si 360degres 

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
def _applyBC(t, metrics, hook1, nstep,omp_mode, var="Density"):
    zones = Internal.getZones(t)
 
    fasts._applyBC(zones, metrics, hook1, nstep, omp_mode, var)
    return None

#==============================================================================
def _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, omp_mode, hook1, nitmax=1, rk=1, exploc=1, num_passage=1): 

   # timecount = numpy.zeros(4, dtype=numpy.float64)

   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None :
           tc_compact = Internal.getNodeFromName1( tc, 'Parameter_real')
           #Si param_real n'existe pas, alors pas de raccord dans tc
           if tc_compact is not  None:

              param_real= tc_compact[1]
              param_int = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1["neq_max"] == 5: varType = 2
              else                    : varType = 21

              for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              type_transfert = 2  # 0= ID uniquement, 1= IBC uniquement, 2= All
              no_transfert   = 1  # dans la list des transfert point a point
              Connector.connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage)#,timecount)

       #apply BC
       #t0=timeit.default_timer()
       if rk != 3 and exploc != 2:
          _applyBC(zones, metrics, hook1, nstep, omp_mode, var=vars[0])
       #t1=timeit.default_timer()
       #print "Time BC",(t1-t0)
            

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
      if var == 'thermique':
        varmy += ['Temperature','T^2','rouT','rovT','rowT','Eps_T' ]
        lgrad =  1

    for i in range(len(varmy)):
      varmy[i] = 'centers:'+varmy[i]

    # on determine le nbr de cellule fictive active pour le calcul des moyennes
    numcellfic = 2 
    ific       = 2   # a adapter en DF
    if lgrad == 1: numcellfic = 1

    zones = []
    for b0 in Internal.getNodesFromType1(t,'CGNSBase_t'):
        if b0[0] != PostBaseName:
            zones += Internal.getZones(b0)

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
            datap[8]  =(dim_my[1]-1)*(dim_my[2]-1)                       # increment k
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

    _compact(tmy, fields=varmy) 

    return tmy

#==============================================================================
# compute statistique in place
#==============================================================================
def _computeStats(t, tmy, metrics):
    """Compute the space/time average of flowfields in tmy."""
    zones    = Internal.getZones(t)
    zones_my = Internal.getZones(tmy)

    bases= Internal.getNodesFromType1(t,'CGNSBase_t')
    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
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

    _compact(tmy, fields=varmy) 

    return tmy

#==============================================================================
# compute enstropy et TKE in place
#==============================================================================
def _computeEnstrophy(t, metrics, time):
    

    zones = Internal.getZones(t)
    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy

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
                 if var=='dDensity/dt'                     : lcompact_drodt = True
                 tmp = numpy.ones( (dim[1]-1,dim[2]-1,dim[3]-1) , dtype=numpy.float64)
                 Internal.createChild(solution, var, 'DataArray_t', tmp)

    if var_loc != []:
       if lcompact_Q:    _compact(zones, fields=['centers:QCriterion'])
       if lcompact_Enst: _compact(zones, fields=['centers:Enstrophy' ])
       if lcompact_Rotx: _compact(zones, fields=['centers:RotX' ])
       if lcompact_drodt: _compact(zones, fields=['centers:dDensitydt' ])

       dtloc = Internal.getNodeFromName3(t , '.Solver#dtloc')  # noeud
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

    var_loc = []
    vargrad = []

    ##verifie si les noeuds existent dans l'arbre
    zones = Internal.getZones(t)
    nd =0
    for z in zones:
        #
        dim   = Internal.getZoneDim(z)
        size = (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        solution = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
     
        for var in vars:
            node = Internal.getNodeFromName1(solution, var)
            #verifie si la variable existe dans l'arbre
            if node is None:
               print("no", var, "in tree: gradient is not computed.")
            else:
               var_loc.append(var)
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

               if lcompact: _compact(z, fields=['centers:'+vargrad[0],'centers:'+vargrad[1],'centers:'+vargrad[2]])
        nd += 1
      
    if var_loc != []:
       dtloc = Internal.getNodeFromName3(t , '.Solver#dtloc')  # noeud
       dtloc = Internal.getValue(dtloc)                       # tab numpy

       # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
       if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)
       fasts.computePT_gradient(zones,  metrics, var_loc, vargrad, FastC.HOOK, order)
    return None

#==============================================================================
# Display
#==============================================================================
def displayTemporalCriteria(t, metrics, nitrun, format=None, gmres=None):
    """Display CFL and convergence information.""" 
    return display_temporal_criteria(t, metrics, nitrun, format, gmres)
    
def display_temporal_criteria(t, metrics, nitrun, format=None, gmres=None):
    zones        = Internal.getZones(t)
    dtloc        = Internal.getNodeFromName3(t, '.Solver#dtloc')
    dtloc_numpy  = Internal.getValue(dtloc)
    nssiter      = int(dtloc_numpy[0])
    nzones	 = len(zones)
 
    #a = Internal.getNodeFromName2(zones[0], 'model')
    #model = Internal.getValue(a)
    #neq = 5
    #if (model == 'nsspalart' or model =='NSTurbulent'): neq = 6
    neq_max = 6
   
    cvg_numpy = numpy.empty((nzones,2*neq_max), dtype=numpy.float64)
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
    residu = fasts.display_ss_iteration( zones, metrics, cvg_numpy, nitrun, nssiter, lft)

    if gmres is None: return None
    else: return residu
#==============================================================================
# Interface for Vtune/Advisor collection control
#==============================================================================
def itt(var):
    if var == 'pause':
          ivar =1
    else :
          ivar = 0
    print("itt collection (Vtune/Advisor)", var)
    fasts.itt(ivar)
    return None
    
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
    

    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy

    zones = Internal.getZones(t)
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it;
    
    bases = Internal.getNodesFromType2(t, 'CGNSBase_t')
    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    fasts.computePT_velocity_ale(zones,  metrics, FastC.HOOK, ompmode)
    return None

#==============================================================================
# init velocity (ALE) from DADS maillage deformable
#==============================================================================
def copy_velocity_ale(t, metrics, it=0):
    

    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy

    zones = Internal.getZones(t)
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it;
    
    bases = Internal.getNodesFromType2(t, 'CGNSBase_t')
    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    #on filtre les zones non deformables
    zones_aleDeformation=[]; metric_aleDeformation=[]
    c=0
    for z in zones:
      tmp = Internal.getNodeFromName1(z,'Motion')
      if tmp is not None:
          zones_aleDeformation.append(z)
          metric_aleDeformation.append(metrics[c])
      c+=1

    fasts.copy_velocity_ale(zones_aleDeformation, metric_aleDeformation , FastC.HOOK, ompmode, it)
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
    varsR   = ['RSD_L2','RSD_oo']
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
        tmp = numpy.empty((nrec), numpy.int32)
        Internal.createChild(c, 'IterationNumber', 'DataArray_t', tmp)
        for var in varsR:
            tmp = numpy.empty((nrec*neq), numpy.float64)
            Internal.createChild(c, var ,'DataArray_t', tmp)
    return None

#==============================================================================
# extraction des residus du Pytree dans fichier tecplot ascii
# ==============================================================================
def extractConvergenceHistory(t, fileout):
    """Extract residuals in an ascii file."""
    zones = Internal.getZonePaths(t)
    it = []; RSD_L2 = []; RSD_oo = []; a = []
    nd = 0
    FileCvgName = fileout
    FileCvg = open(FileCvgName,'w')
    for i in zones:
        node = Internal.getNodeFromPath(t, i+'/ZoneConvergenceHistory')
        lastRec = node[1][0]
        #it = C.convertPyTree2Array(i+"/ZoneConvergenceHistory/IterationNumber", t)
        #a  = C.convertPyTree2Array(i+"/ZoneConvergenceHistory/RSD_L2", t)
        it = Internal.getNodeFromPath(t,i+"/ZoneConvergenceHistory/IterationNumber")
        a = Internal.getNodeFromPath(t, i+"/ZoneConvergenceHistory/RSD_L2")
        nrec = it[1].size
        neq = a[1].size//nrec
        RSD_L2 = numpy.reshape(a[1],(neq,nrec),order='F')
        #a = C.convertPyTree2Array(i+"/ZoneConvergenceHistory/RSD_oo", t)
        a = Internal.getNodeFromPath(t, i+"/ZoneConvergenceHistory/RSD_oo")
        RSD_oo = numpy.reshape(a[1],(neq,nrec),order='F')
        a='"'
        if neq == 5: 
            var="VARIABLES = it RO_l2 ROU_l2 ROV_l2 ROW_l2 ROE_l2 ROoo ROUoo ROVoo ROWoo ROEoo\n"
        if neq == 6: 
            var="VARIABLES = it RO_l2 ROU_l2 ROV_l2 ROW_l2 ROE_l2 NUT_l2 ROoo ROUoo ROVoo ROWoo ROEoo NUToo\n"
        if nd == 0: FileCvg.write("%s"%(var)) 
        nd =+ 1
        FileCvg.write("ZONE T=%sbloc %s %s I=%d F=POINT\n"%(a,i,a,lastRec))
        c = it[1]
        for l in range(lastRec):
           a  = ""
           for k in range(neq): a = a+"{0:7f} ".format(RSD_L2[(k,l)])
           for k in range(neq): a = a+"{0:7f} ".format(RSD_oo[(k,l)])
           FileCvg.write('%s %s\n'%(1+c[l],a))
    FileCvg.close()

#==============================================================================
# Cree un arbre Stress pour calcul effort
# IN: t: tree
# IN: BC: list des BC concernees ['BCWall', BCFarfield',...]
# IN: window: ou range de window
# OUT: return arbre stress
#==============================================================================
def createStressNodes(t, BC=None, windows=None):
    """Create nodes to store stress data."""
    import Converter.GhostCells as Ghost
    try: import Transform.PyTree as T
    except: raise ImportError("createStressNodes: requires transform module.")

    PostBaseName = 'Stress' # nom de la base POST
    vars0  = ['CoordinateX','CoordinateY','CoordinateZ']
    var    = ['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity']

    teff = C.newPyTree([PostBaseName])
    b    = Internal.getNodesFromName1(teff, PostBaseName)

    zones = []
    no_z = 0
    for b0 in Internal.getNodesFromType1(t,'CGNSBase_t'):
        dimbase=  b0[1][0]
        if b0[0] != PostBaseName:
           zones = Internal.getNodesFromType1(b0, 'Zone_t')

           for z in zones:

              # gestion zone2d
              dimzone= Internal.getZoneDim(z)
              inck0 = 2
              inck1 = 1
              vargrad=['gradxVelocityX','gradyVelocityX','gradzVelocityX','gradxVelocityY','gradyVelocityY','gradzVelocityY','gradxVelocityZ','gradyVelocityZ','gradzVelocityZ','gradxTemperature','gradyTemperature','gradzTemperature', 'CoefPressure','ViscosityMolecular']
              if dimzone[3] == 2: 
                   inck0=0
                   inck1=0
                   vargrad=['gradxVelocityX','gradyVelocityX','gradxVelocityY','gradyVelocityY','gradxTemperature','gradyTemperature','CoefPressure','ViscosityMolecular']

              varc  = []
              for i in range(len(var)): varc.append('centers:'+var[i])
              for i in range(len(vargrad)): varc.append('centers:'+vargrad[i])
              
              if BC is not None:
                 bc = Internal.getNodeFromName1(z, "ZoneBC")
                 list_bc =[]
                 if bc is not None: list_bc = bc[2]
              else:
                 list_bc =['window']
              ific    = 2
              param_int = Internal.getNodeFromName2(z, 'Parameter_int')
              if param_int is not None: ific = param_int[1][3]
              ndf =0
              for v in list_bc:
                 if BC is not None:
                   name = Internal.getValue(v)
                 else:
                   name = v[0]

                 if windows is None:
                     windows = [None]

                 inum = 0

                 for window in windows:

                   if (BC is not None and name in BC) or (window is not None and z[0]==window[0]):
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
                        ptrg = [ 'PointRange', wrange,  [], 'IndexRange_t']
                      
                     dim = Internal.getValue(ptrg)
                     
                     #dir=0,...5
                     idir = Ghost.getDirection__(dimbase, [ptrg])
                     if windows[0] is not None:
                       if   window[7]== 'imin': idir =0
                       elif window[7]== 'imax': idir =1
                       elif window[7]== 'jmin': idir =2
                       elif window[7]== 'jmax': idir =3
                       elif window[7]== 'kmin': idir =4
                       elif window[7]== 'kmax': idir =5
                     
                     inci =0
                     incj =0
                     inck =0
                     if BC is not None:
                       if  (idir==0): inci= ific
                       elif(idir==1): inci=-ific
                       elif(idir==2): incj= ific
                       elif(idir==3): incj=-ific
                       elif(idir==4): inck= inck0
                       elif(idir==5): inck=-inck0
  
                     ni = dim[0,1]-dim[0,0]
                     nj = dim[1,1]-dim[1,0]
                     nk = dim[2,1]-dim[2,0]

                     ideb = dim[0,0]+inci
                     ifin = dim[0,1]+inci
                     jdeb = dim[1,0]+incj
                     jfin = dim[1,1]+incj
                     kdeb = dim[2,0]+inck
                     kfin = dim[2,1]+inck

                     #print 'subZ', ideb,ifin,jdeb,jfin,kdeb,kfin

                     zp = T.subzone(z, (ideb,jdeb,kdeb), (ifin,jfin,kfin))
                     zp[0] = z[0]+'_'+v[0]+str(inum)
                     inum  = inum + 1
                     C._extractVars(zp, vars0)
                     c = 0
                     for v1 in varc: 
                        C._initVars(zp, v1, c)
                        c += 1
                     #print 'zone name=',zp[0]
                     #print 'dim0=',dim[0,0],dim[0,1],dim[1,0],dim[1,1],dim[2,0],dim[2,1]

                     if idir <= 1:
                       ni1  = nj
                       nj1  = nk
                       jmin = jdeb -ific      #jmin
                       jmax = jfin -ific-1    #jmax
                       kmin = kdeb -inck0     #kmin
                       kmax = kfin -inck0-1   #kmax
                       imin = ideb -ific
                       imax = imin
                     elif(idir <= 3):   
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
                        raise ValueError("_createStatNodes: Parameter_int is missing for zone %s."%zp[0])
 
                     #print 'loop effort',imin,imax,jmin,jmax,kmin,kmax
                     ##modifie le moeud param_int de l'arbre teffot (issue de subzone) pour la fonction initNuma
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
                        
                     #print 'dim1=',imin,imax,jmin,jmax,kmin,kmax
                     #print 'idir=',idir+1

                     _compact(zp, fields=varc, init=False) #allocation compact uniquememnt; init dans er calcul effort

                     b[0][2].append(zp)


              no_z +=1
              ndf  +=1

    Internal._rmNodesByType(b,'ZoneGridConnectivity_t')
    Internal._rmNodesByType(b,'ZoneBC_t')
    Internal._rmNodesByType(b,'Rind_t')
    Internal._rmNodesByName(b,'.Solver#define')
    Internal._rmNodesByName(b,'Parameter_real')
    Internal._rmNodesByName(b,'CFL_minmaxmoy')
    Internal._rmNodesByName(b,'type_zone')
    Internal._rmNodesByName(b,'model')
    Internal._rmNodesByName(b,'temporal_scheme')
    Internal._rmNodesByName(b,'GridCoordinates#Init')
    Internal._rmNodesByName(b,'FlowSolution')

    return teff

#==============================================================================
#
# Calcul des effort (in place)
#
#==============================================================================
def _computeStress(t, teff, metrics, xyz_ref=(0.,0.,0.) ):
    """Compute efforts in teff."""
    
    zones     = Internal.getZones(t)
    zones_eff = Internal.getZones(teff)
    
    bases= Internal.getNodesFromType1(t,'CGNSBase_t')
    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)
    
    # Cree des tableaux temporaires de travail (wiggle, coe, drodm, lok, iskip_lu)
    if FastC.HOOK is None: 
            dtloc  = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
            dtloc  = Internal.getValue(dtloc)                       # tab numpy
            FastC.HOOK   = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT)
            nitrun =0; nstep =1

            distrib_omp = 0
            hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )

    else:  hook1  = FastC.HOOK

    effort  = numpy.empty(8, numpy.float64)
    pos_eff = numpy.empty(3, numpy.float64)
    pos_eff[0] = xyz_ref[0]; pos_eff[1] = xyz_ref[1]; pos_eff[2] = xyz_ref[2]

    #print 'ompmode=', ompmode
    fasts.compute_effort(zones, zones_eff, metrics, hook1, effort, pos_eff, ompmode)

    return effort

#==============================================================================
#
# Calcul distributiona threads
#
#==============================================================================
def distributeThreads(t, metrics, work, nstep, nssiter, nitrun, Display=False):

  zones          = Internal.getZones(t)
  mx_omp_size_int= work["MX_OMP_SIZE_INT"]
  #print 'mx_omp_size_int', mx_omp_size_int
  display = 0
  if Display: display = 1

  fasts.distributeThreads(zones , metrics, nstep, nssiter , nitrun, mx_omp_size_int, display)

  return None


#==============================================================================
#
# IBC ordre 0: preparation
# surf: arbre avec une base par obstacle (plusieurs zones possibles)
#==============================================================================
def setIBCData_zero(t, surf, dim=None):

  zones = Internal.getZones(t)

  # recup des obstacles
  bodies = []
  for b in Internal.getBases(surf):
      body = Internal.getZones(b)
      bodies.append(body)
  
  # Matrice de masquage (arbre d'assemblage)
  nbodies = len(bodies)
  BM = numpy.ones((1,nbodies),dtype=numpy.int32)
  if dim is None:
     sol= Internal.getNodeFromName1(zones[0],'FlowSolution#Centers')
     nk = Internal.getNodeFromName1(sol,'Density')[1].shape[2]
     if nk ==1 : dim = 2
     else:       dim = 3

  C._initVars(t,"{centers:cellN_IBC}=1.")
  #X._blankCells(t, bodies, BM, blankingType="center_in", dim= dim,delta=0., cellNName='cellN_IBC')
  t = X.blankCells(t, bodies, BM, blankingType="center_in", dim= dim,delta=0., cellNName='cellN_IBC')

  numz = {}
  zones = Internal.getZones(t)
  for z in zones:
     cellN_node = Internal.getNodeFromName2(z,'cellN_IBC')
     #cellN_node[0]= 'cellN_IBC'
     cellN_IBC = cellN_node[1]

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

        print('zone', z[0],': plage IBC=', ibc[1:7])
      
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

   z[0]=zname
   zout.append(z)

 tp = C.newPyTree([basename])
 # Attache les zones
 tp[2][1][2] += zout
  
 return tp


#==============================================================================
# display CPU efficiency diagnostic
#==============================================================================
def display_cpu_efficiency(t, mask_cpu=0.08, mask_cell=0.01, diag='compact', FILEOUT='listZonesSlow.dat', RECORD=None):

 bases = Internal.getNodesFromType1(t      , 'CGNSBase_t')       # noeud
 own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
 dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

 node = Internal.getNodeFromName1(bases[0], '.Solver#define')
 node = Internal.getNodeFromName1(node, 'omp_mode')
 ompmode = OMP_MODE
 if  node is not None: ompmode = Internal.getValue(node)

 dtloc = Internal.getValue(dtloc) # tab numpy
 ss_iteration  = int(dtloc[0]) 

 timer_omp = FastC.HOOK["TIMER_OMP"]
 ADR = OMP_NUM_THREADS*2*(ss_iteration)
 echant    =  timer_omp[ ADR ]
 if echant == 0.:
   print('nombre iterations insuffisant pour diagnostic: nitrun * ss_iteration > 15')
   return None

 zones = Internal.getZones(t)
 tps_percell     =0.
 tps_percell_max =0.
 cout            =0.
 cells_tot       =0
 data            ={}
 for z in zones:
    echant    =  timer_omp[ ADR ]
    param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
    solver_def= Internal.getNodeFromName2(z, '.Solver#define')
    if ompmode ==1:
       Ptomp       = param_int[69]
       PtrIterOmp  = param_int[Ptomp] 
       PtZoneomp   = param_int[PtrIterOmp]
       NbreThreads = param_int[ PtZoneomp  + OMP_NUM_THREADS ]
    else:
       NbreThreads = OMP_NUM_THREADS

    if diag== 'compact':
      tps_zone_percell     =0.
      tps_zone_percell_max =0.
      ithread_max          = 0
      for i in range(OMP_NUM_THREADS):
          if ompmode ==1:
            ithread = param_int[ PtZoneomp  +  i ]
            if ithread != -2:
               tps_zone_percell += timer_omp[ ADR + 1+ithread ]
               if timer_omp[ ADR + 1+ithread ] > tps_zone_percell_max: 
                    tps_zone_percell_max = timer_omp[ ADR + 1+ithread ]
                    ithread_max          = ithread
          else:
            tps_zone_percell += timer_omp[ ADR + 1+i ]
            if timer_omp[ ADR + 1+i ] > tps_zone_percell_max: 
                 tps_zone_percell_max = timer_omp[ ADR + 1+i ]
                 ithread_max          = i
      
      ijkv             = param_int[20]*param_int[21]*param_int[22]
      tps_percell     += tps_zone_percell/echant/NbreThreads*ijkv
      tps_percell_max += tps_zone_percell_max/echant*ijkv
      cout_zone        = tps_zone_percell/echant/NbreThreads*ijkv
      cout            += cout_zone
      data[ z[0] ]   = [ tps_zone_percell/echant/NbreThreads, cout_zone]

      cells_tot   += ijkv
      if RECORD is None: print('zone= ', z[0],'cpumoy/cell= ',tps_zone_percell/echant/NbreThreads,'cpumax= ',tps_zone_percell_max/echant,'thread lent=', ithread_max ,'Nthtread actif=', NbreThreads, ' dim zone=',ijkv,'dim tot=',cells_tot)
      if RECORD is not None:
          tape = tps_zone_percell/echant/NbreThreads

      tps_zone_percell = max(tps_zone_percell, 1.e-11)
      tps_zone_percell_max = max(tps_zone_percell_max, 1.e-11)
      
      perfo = numpy.empty(2, dtype=numpy.float64)
      perfo[0]= int(echant*NbreThreads/tps_zone_percell)
      perfo[1]= int(echant/tps_zone_percell_max);
      Internal.createUniqueChild(solver_def, 'Cups', 'DataArray_t', value=perfo)

    else:
      for i in range(OMP_NUM_THREADS):
          if ompmode ==1:
            ithread = param_int[ PtZoneomp  +  i ]
            if ithread != -2: 
               if RECORD is None: print('zone= ', z[0], 'cpu= ',timer_omp[ ADR + 1+ithread ]/echant,' th=  ', ithread, 'echant= ', echant)
          else:
            ithread = i
            if RECORD is None: print('zone= ', z[0], 'cpu= ',timer_omp[ ADR + 1+ithread ]/echant,' th=  ', ithread, 'echant= ', echant)

    ADR+= OMP_NUM_THREADS+1

 tps_percell/=cells_tot
 if RECORD is None: print('cpu moyen %cell en microsec: ', tps_percell, tps_percell_max/cells_tot)

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

#==============================================================================
# compute rk3 dtloc 
# 
#==============================================================================
def _computeguillaume1(t, metrics, nitrun, tc=None, graph=None, layer="Python", NIT=1):

    

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t') 

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    exp_local = Internal.getNodeFromName1(node, 'exp_local') 
    explocal = Internal.getValue(exp_local)
    r_k = Internal.getNodeFromName1(node, 'rk') 
    rk = Internal.getValue(r_k)

    
    dtloc = Internal.getValue(dtloc) # tab numpy
    nitmax = int(dtloc[0])   
   

    #### a blinder...
    itypcp = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1][29]
    #### a blinder...


    if tc is not None:
         bases = Internal.getNodesFromType1(tc, 'CGNSBase_t')  # noeud
         tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
         #print tc_compact
         if tc_compact is not None:
                param_real= tc_compact[1]
                param_int = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]
       
                zonesD = []
                for f in bases:
                    tmp = Internal.getNodesFromType1(f, 'Zone_t') 
                    zonesD += tmp

    zonesD    = Internal.getZones(tc)

    nbcomIBC    = param_int[1]
    shift_graph = nbcomIBC + param_int[2+nbcomIBC] + 2
    comm_P2P = param_int[0]
    pt_ech = param_int[comm_P2P + shift_graph]
    dest   = param_int[pt_ech]


    #nitmaxs=8
    if nitrun == 1: print('Info: using layer trans=%s (ompmode=%d)'%(layer, ompmode))

    if layer == 'Python':
           
        for nstep in range(1, nitmax+1): # Etape explicit local

            # determination taille des zones a integrer (implicit ou explicit local)

            hook1  = FastC.HOOK.copy()
            distrib_omp = 0
            hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )
            nidom_loc = hook1["nidom_tot"]
            
            skip      = 0
            if (hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1): skip = 1


            # calcul Navier_stokes + appli CL
            if nidom_loc > 0 and skip == 0:

                # Navier-Stokes
                nstep_deb = nstep
                nstep_fin = nstep
                layer_mode= 0
                nit_c     = 1
                fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, hook1)

                fasts.dtlocal2para_(zones,zonesD,param_int,param_real,hook1,0,nstep,ompmode,1,dest)
                
                if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
                elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
                _applyBC(zones,metrics, hook1, nstep, ompmode, var=vars[0])    

                FastC.switchPointers2__(zones,nitmax,nstep)
                
                # Ghostcell
                if (nitmax%3 != 0) : # Tous les schemas sauf constantinescu RK3
                    if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
                    elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1']
                    timelevel_target = int(dtloc[4])
                    _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, nitmax, rk, explocal)
            
                fasts.recup3para_(zones,zonesD,param_int,param_real,hook1,0,nstep,ompmode,1) 

                if (nstep%2==0) :
                    timelevel_target = int(dtloc[4])
                    vars = ['Density'  ]
                    _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, nitmax, hook1, nitmax, rk, explocal, 2)

                if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
                elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
                _applyBC_(zones,metrics, hook1, nstep, ompmode,  var=vars[0])    

    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT
      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c , FastC.HOOK)

      FastC.switchPointers3__(zones,nitmax)
                         
    # switch pointers
    if (layer == 'Python') :
        FastC.switchPointers__(zones, 1, 3)

    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT = 1


    return None


#==============================================================================
# compute guillaume2
#==============================================================================
def _computeguillaume2(t, metrics, nitrun, tc=None, graph=None, layer="Python", NIT=1):

    

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t') 




    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    exp_local = Internal.getNodeFromName1(node, 'exp_local') 
    explocal = Internal.getValue(exp_local)
    r_k = Internal.getNodeFromName1(node, 'rk') 
    rk = Internal.getValue(r_k)


    dtloc = Internal.getValue(dtloc) # tab numpy
    nitmax = int(dtloc[0])   
    #print('nitmax= %d'%nitmax)
 
    #### a blinder...
    itypcp = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1][29]
    #### a blinder...

    #param_intcopy=[]
    #param_realcopy=[]
    if tc is not None:
         bases = Internal.getNodesFromType1(tc, 'CGNSBase_t')  # noeud
         tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
         if tc_compact is not None:
                param_real= tc_compact[1]
                param_int = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]
 
         
                zonesD = []
                for f in bases:
                    tmp = Internal.getNodesFromType1(f, 'Zone_t') 
                    zonesD += tmp

    taille_tabs = 20000000
    rostk = numpy.empty(taille_tabs, dtype=numpy.float64)  # tableau de stockage des valeurs 
    drodmstk = numpy.empty(taille_tabs, dtype=numpy.float64) # tableau de stockage des flux (drodm)
    constk = numpy.empty(taille_tabs, dtype=numpy.float64) # tableau de stockage des flux (conservativite)

    zonesD    = Internal.getZones(tc)



    #print('nitmax= %d'%nitmax)
    #nitmaxs=1

    if nitrun == 1: print('Info: using layer trans=%s (ompmode=%d)'%(layer, ompmode))

    if layer == "Python": 
    
        for nstep in range(1, nitmax+1): # Etape explicit local
 
            #print('nstep '%nstep)             
            # determination taille des zones a integrer (implicit ou explicit local)

            hook1  = FastC.HOOK.copy()
            distrib_omp = 0
            hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )
            nidom_loc = hook1["nidom_tot"]

            skip      = 0
            if (hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1): skip = 1

            # calcul Navier_stokes + appli CL
            if nidom_loc > 0 and skip == 0:

                # Navier-Stokes
                nstep_deb = nstep
                nstep_fin = nstep
                layer_mode= 0
                nit_c     = 1
                fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, hook1) 

                #apply BC apres la fonction recup
                if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
                elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
                _applyBC(zones,metrics, hook1, nstep, ompmode, var=vars[0])
           
                #fonction qui controle stockage, destockage et calcule le predicteur sans passer par le code entier
                fasts.stockrecup(zones,rostk,drodmstk,constk,hook1,nstep,taille_tabs)

                # Ghostcell
                if (nitmax%3 != 0) : # Tous les schemas sauf constantinescu RK3
                    if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
                    elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1']
                    timelevel_target = int(dtloc[4])
                    _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1) 

                elif (nitmax%3 == 0) : # Constantinescu RK3             
                    if   (nstep%3==1)  and itypcp == 2 :  vars = ['Density_P1']
                    elif (nstep%3==2)  and itypcp == 2 :  vars = ['Density_M1']
                    elif (nstep%3==0)  and itypcp == 2 :  vars = ['Density'  ]
                    timelevel_target = int(dtloc[4])
                    _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1)
    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT
      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c , FastC.HOOK)

        
    # switch pointers
    FastC.switchPointers__(zones,1,2)
    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT = 1

    return None


#==============================================================================
# temps d'estimation rk
#==============================================================================
def _temps_estim(t,varN,varP,drodm_,coe_,nstep):
    #if (nstep == 1):
    zones = Internal.getZones(t)
    #else:
    #    tp = C.convertFile2PyTree('lamb__.cgns')
    #    zones = Internal.getZones(tp)        
        
    index = numpy.empty(len(zones)+1, dtype=numpy.int32)
    index_ = numpy.empty(len(zones)+1, dtype=numpy.int32)
    ind=0
    ind_=0
    i=0
    for z in zones:
        dim = Internal.getZoneDim(z)
        ind = ind + 5*(dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        ind_= ind_ + (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        index[i+1]=ind
        index_[i+1]=ind_
        i=i+1
    index[0]=0
    index_[0]=0

    # Nouvel arbre
    tp = C.newPyTree(['Base'])

    i=0
    out=[]
    for z in zones:
         #if (nstep == 1):
         z = C.addVars(z, 'centers:temps')
         z = C.initVars(z, 'centers:temps', 0.)
         dim = Internal.getZoneDim(z)
         node   = Internal.getNodeFromName1(z ,'FlowSolution#Centers')
         var_P1 = Internal.getNodeFromName1(node, varP[0])
         var_N  = Internal.getNodeFromName1(node, varN[0])
         ndimdx = (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
         drodm  = numpy.copy(drodm_[index[i]:index[i] + ndimdx])
         coe    = numpy.copy(coe_[index_[i]:index_[i] + ndimdx])
         #drodm = numpy.reshape(drodm,(dim[1]-1,dim[2]-1,dim[3]-1))
         for j in range(2,dim[2]-3):
            for k in range(2,dim[1]-3):
                 diff =  abs(var_N[1][k,j,0] - var_P1[1][k,j,0])/drodm[k + j*(dim[1]-1)]
                 #diff = abs(diff)/coe[k + j*(dim[1]-1)]
                 #a = C.getValue(z,'centers:temps', (k,j,0))
                 C.setValue(z, 'centers:temps', (k,j,0),diff)
         i = i + 1
         out.append(z)

    # Attache les zones
    tp[2][1][2] += out

    return tp 
#==============================================================================
# Calcule la CFL en chaque maille et la met dans l'arbre t
#==============================================================================
def _comp_cfl(t,metrics):

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t') 
    zones = []
    for b in bases:
        zones += Internal.getNodesFromType1(b, 'Zone_t')

    omp_mode = 0

    fasts.prep_cfl(zones, metrics, 1, 1, omp_mode)
    return None
#==============================================================================
# decoupe maillage pour dtloc instationnaire
#==============================================================================
def _decoupe(t):

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t') 
    zones = []
    for b in bases:
        zones += Internal.getNodesFromType1(b, 'Zone_t')

    infos = fasts.decoupe_maillage(zones,4)
    return infos
#==============================================================================
# compute_dpJ_dpW in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute_dpJ_dpW(t, teff, metrics, cosAoA, sinAoA, surfinv):

    

    zones     = Internal.getZones(t)
    zones_eff = Internal.getZones(teff)

    if FastC.HOOK is None: 
            dtloc  = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
            dtloc  = Internal.getValue(dtloc)                       # tab numpy
            FastC.HOOK   = FastC.createWorkArrays__(zones, dtloc, FastC.FIRST_IT); 
            nitrun =0; nstep =1;
            hook1  = FastC.HOOK.copy()
            distrib_omp = 0
            hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )
    else:   hook1  = FastC.HOOK

    fasts.compute_dpJ_dpW(zones, zones_eff, metrics, hook1, cosAoA, sinAoA, surfinv)
    return None


#==============================================================================
# computeAdjoint in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _computeAdjoint(t, metrics, nit_adjoint, indFunc, tc=None, graph=None):

    

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    zones = []
    for f in bases:
        zones += Internal.getNodesFromType1(f, 'Zone_t') 

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = OMP_MODE
    if  node is not None: ompmode = Internal.getValue(node)

    dtloc = Internal.getValue(dtloc) # tab numpy
    nitmax = int(dtloc[0])                 
    orderRk = int(dtloc[len(dtloc)-1])
    

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

#    nstep=1
#    hook1  = FastC.HOOK.copy()
#    hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nit_adjoint, nstep) )
#    hook1  = FastC.HOOK
    

   #--------------------------------------------------------------------------
   #  0  peut-on faire ici un test pour verifier que dpJdpW a bien ete prealablement
   #   calcule (ordre correct des appelants) (boolean isdpJpWComputed ... ?) IVAN
   # champ adjoint et increment adjoint initialise a 0 (?) a verifier
   #--------------------------------------------------------------------------------

  # if nit_adjoint == 1 :
       # if indFunc == 1 :
       #      vars = ['dpCLp_dpDensity']
       #      Connector.connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, varType, type_transfert, no_transfert)

       # if indFunc == 2 : 
       #      vars = ['dpCDp_dpDensity']
       #      Connector.connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, varType,  type_transfert, no_transfert)

       # if vars is not None:
            # compute_dpJ_dpW = True       #dpJ_dpW is already computed
            # print "dpJ_dpW is already computed"
       # else:
            # print "dpJ_dpW is missing"   #dpJ_dpW is not computed
            # return None                  #This is equivalent to exit(), because dpJ_dpW has not been computed, so we can not compute RhsIter_Adjoint.

   #--------------------------------------------------------------------------------
   #  1 raccord sur l'adjoint via les ghost cells
   #-----------------------------------

#    if indFunc == 1:
#        vars = ['AdjCLp_RDensity']     # should work for the five variables of incAdj (?)
#           
#    if indFunc == 2: 
#        vars = ['AdjCDp_RDensity']     # should work for the five variables of incAdj (?)
#
#    # apply transfers
#    if tc is not None and hook1[12] ==0: 
#       if   hook1[10] == 5: varType = 2
#       else               : varType = 21

#    #print 'transfert', nstep, skip,hook1[13], hook1[12]        
#    if tc_compact is not None:
#       for v in vars: C._cpVars(t, 'centers:'+v, tc, v)
#       type_transfert = 2  # 0= ID uniquement, 1= IBC uniquement, 2= All
#       no_transfert   = 1  # dans la list des transfert point a point

#    Connector.connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, varType, type_transfert, no_transfert)


#    #-----------------------------------------------------------------
#    #  2 calcul membres de droite et de gauche de l'algo iteratif pour l'adjoint 
#    #--------------------------------------------

#    fasts.compute_RhsIterAdjoint(zones, metrics, nit_adjoint, nstep, indFunc, omp_mode, hook1)

#    return None
 

    #-----------------------------------------------------------------
    #  3 calcul membres de droite de l'algo iteratif pour l'adjoint 
    #--------------------------------------------
    
''' 
     vars = ['IncAdj']     

     do i=1,6 
  
       fasts.compute_LorUAdjoint(zones, metrics, nitrun, nstep, indFunc, omp_mode, hook1)

       Connector.connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, varType, type_transfert, no_transfert)
       enddo 
    #-----------------------------------------------------------------
    #  4 update
    #------------------------------------------------------
 
    if indFunc == 1:
        # adjCLp = adjCLp + incAdj
    if indFunc == 2: 
        # adjCDp = adjCDp + incAdj
'''

   # return None

'''
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
    orderRk = int(dtloc[len(dtloc)-1])
    

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

