# FastS + MPI
import PyTree
import fasts
from PyTree import display_temporal_criteria, createConvergenceHistory, extractConvergenceHistory, createStressNodes, _computeStress, createPrimVars, _createPrimVars, createStatNodes, _computeStats, initStats, _computeEnstrophy, _computeVariables, _computeGrad, _compact, _applyBC, _buildOwnData, _init_metric, allocate_metric, _BCcompact, _motionlaw, _movegrid, _computeVelocityAle, checkBalance, itt, HOOK, distributeThreads, _build_omp
import timeit

import numpy
import sys

try:
    import Converter.PyTree as C
    import Converter.Mpi    as Cmpi
    import Distributor2.PyTree as D2
    import Converter.Internal as Internal
    import Connector.Mpi as Xmpi
    import Connector.PyTree as X
    import Connector
    import Fast.Internal as FastI
    import os
    import math
except:
    raise ImportError("FastS: requires Converter and Connector modules.")

#==============================================================================
#def _compute(t, metrics, nitrun, tc=None, graph=None, layer="Python"):
def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c"):
    if graph is not None:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
    else: 
        procDict=None; graphID=None; graphIBCD=None

    base = Internal.getNodeFromType1(t,"CGNSBase_t")
    own   = Internal.getNodeFromName1(base, '.Solver#ownData')  
    dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')
    zones = Internal.getZones(t)
    node = Internal.getNodeFromName(t, '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    omp_mode = 0
    if node is not None: omp_mode = Internal.getValue(node)
    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])
    orderRk = int(dtloc[len(dtloc)-1])

    #### a blinder...
    itypcp = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1][29]
    #### a blinder...

    if layer == "Python": 

      for nstep in xrange(1, nitmax+1): # pas RK ou ssiterations

         hook1 = PyTree.HOOK.copy()
         hook1.update(  fasts.souszones_list(zones, metrics, PyTree.HOOK, nitrun, nstep) )
         nidom_loc = hook1["nidom_tot"]

         skip = 0
         if (hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1): skip = 1

         # calcul Navier Stokes + appli CL
         if nidom_loc > 0 and skip ==0:
            # Navier-Stokes
            nstep_deb = nstep
            nstep_fin = nstep
            layer_mode= 0
            fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, omp_mode, hook1)

            # Ghostcell
            vars = PyTree.varsP
            if  nstep == 2 and itypcp == 2 : vars = PyTree.varsN  # Choix du tableau pour application transfer et BC
            timelevel_target = int(dtloc[4])
            _fillGhostcells(zones, tc, metrics, timelevel_target , vars, nstep, omp_mode, hook1, graphID, graphIBCD, procDict)

    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, omp_mode, PyTree.HOOK)

    # data update for unsteady joins
    dtloc[3] +=1   #time_level_motion
    dtloc[4] +=1   #time_level_target

    # switch pointers
    FastI.switchPointers__(zones, orderRk)
    # flag pour derivee temporelle 1er pas de temps implicit
    PyTree.HOOK["FIRST_IT"]  = 1
    PyTree.FIRST_IT          = 1
    return None
#==============================================================================
def _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, omp_mode, hook1, graphID, graphIBCD, procDict): 

   rank=Cmpi.rank 

   #timecount = numpy.zeros(4, dtype=numpy.float64)
   timecount = []
   
   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None:
           tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
           if tc_compact is not None:
              param_real = tc_compact[1]
                
              param_int = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1['neq_max'] == 5: varType = 2
              else                    : varType = 21

              bcType = PyTree.HOOKIBC[0]; Gamma=PyTree.HOOKIBC[1]; Cv=PyTree.HOOKIBC[2]; Mus=PyTree.HOOKIBC[3]; Cs=PyTree.HOOKIBC[4]; Ts=PyTree.HOOKIBC[5]

              if nstep <= 3: 
                 for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              # #recuperation Nb pas instationnnaire dans tc
              type_transfert = 1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        bcType=bcType, varType=varType, compact=1,
                                        graph=graphIBCD, procDict=procDict, 
                                        Gamma=Gamma, Cv=Cv, MuS=Mus, Cs=Cs, Ts=Ts)
              type_transfert = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        varType=varType, compact=1, graph=graphID, procDict=procDict)
       # if (rank == 0):
       #     print "Time in MPI send_buffer, irecv ","%.6f"%timecount[0]
       #     print "Time InterpTransfert (Inter)  ","%.6f"%timecount[1]
       #     print "Time InterpTransfert (Intra)  ","%.6f"%timecount[2]
       #     print "Time in getTransfersInter ","%.6f"%timecount[3]
       # if (rank == 0 ): t0=timeit.default_timer()
       #apply BC
       _applyBC(zones, metrics, hook1, nstep, omp_mode, var=vars[0])
       # if (rank == 0 ):
       #     t1=timeit.default_timer()
       #     print "time/it (s) BC only (=t_RK X 1.0): ",(t1-t0)

   return None
#==============================================================================
def warmup(t, tc, graph=None, infos_ale=None, Adjoint=False, tmy=None):
    #global PyTree.FIRST_IT, PyTree.HOOK, PyTree.HOOKIBC

    rank = Cmpi.rank

    if graph is not None:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
    else: 
        procDict=None; graphID=None; graphIBCD=None

    # Get omp_mode
    omp_mode = 0
    node = Internal.getNodeFromName2(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: omp_mode = Internal.getValue(node)

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastI._reorder(t, tc, omp_mode)

    # Construction param_int et param_real des zones
    _buildOwnData(t)

    # determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: xrange(22,1000)
    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy
    zones = Internal.getZones(t)
    f_it = PyTree.FIRST_IT
    if PyTree.HOOK is None: PyTree.HOOK = FastI.createWorkArrays__(zones, dtloc, f_it ); PyTree.FIRST_IT = f_it

    # allocation d espace dans param_int pour stockage info openmp
    _build_omp(t) 

    # alloue metric: tijk, ventijk, ssiter_loc
    # init         : ssiter_loc
    metrics = allocate_metric(t)

    # Contruction BC_int et BC_real pour CL
    _BCcompact(t) 

    #determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: xrange(22,1000)
    for nstep in xrange(1, int(dtloc[0])+1):
        hook1 = PyTree.HOOK.copy()
        hook1.update(  fasts.souszones_list(zones, metrics, PyTree.HOOK, 1, nstep) )
        distributeThreads(t, metrics, PyTree.HOOK, nstep, int(dtloc[0]) )

    _init_metric(t, metrics, hook1, omp_mode)

    # compact + align + init numa
    rmConsVars=True
    adjoint   =Adjoint
    t = createPrimVars(t, omp_mode, rmConsVars, adjoint)

    zones = Internal.getZones(t) # car create primvar rend zones caduc

    #Allocation HOOKIBC
    if PyTree.HOOKIBC is None: PyTree.HOOKIBC = FastI.getIBCInfo__(t)

    if PyTree.HOOK["neq_max"] == 5: varType = 2
    else                          : varType = 21

    PyTree.HOOK['param_int_ibc'][0] = varType
    PyTree.HOOK['param_int_ibc'][1] = PyTree.HOOKIBC[0]
    PyTree.HOOK['param_real_ibc'][0]= PyTree.HOOKIBC[1]
    PyTree.HOOK['param_real_ibc'][1]= PyTree.HOOKIBC[2]
    PyTree.HOOK['param_real_ibc'][2]= PyTree.HOOKIBC[3]
    PyTree.HOOK['param_real_ibc'][3]= PyTree.HOOKIBC[4]
    PyTree.HOOK['param_real_ibc'][4]= PyTree.HOOKIBC[5]

    #corection pointeur ventijk si ale=0: pointeur Ro perdu par compact.
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
        else: ale =True
        c += 1

    #
    # mise a jour vitesse entrainememnt
    #
    if ale == True and infos_ale is not None:
        print "ale actif. Teta et tetap=", infos_ale
        teta = infos_ale[0];  tetap = infos_ale[1]
        _motionlaw(t, teta, tetap)
        _computeVelocityAle(t,metrics)
    #
    # Compactage arbre transfert
    #
    if tc is not None:      
       X.miseAPlatDonorTree__(zones, tc, graph=graph)

       PyTree.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int')[1]
       param_real_tc               = Internal.getNodeFromName1( tc, 'Parameter_real')
       if param_real_tc is not None: PyTree.HOOK['param_real_tc']= param_real_tc[1]
    else:
        PyTree.HOOK['param_real_tc'] = None
        PyTree.HOOK['param_int_tc']  = None 

    #
    # Compactage arbre moyennes stat
    #
    if tmy is not None:
        sol = Internal.getNodesFromName3(tmy, 'FlowSolution#Centers')
        var = Internal.getNodesFromType1(sol[0] , 'DataArray_t')
        varmy=[]
        for v in var: varmy.append('centers:'+v[0])
        _compact(tmy, fields=varmy)

    #
    # remplissage ghostcells
    #
    hook1['skip_lu'] = 0
    nstep            = 1
    nitrun           = 0
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4]) 
    _fillGhostcells(zones, tc, metrics, timelevel_target, ['Density'], nstep, omp_mode, hook1,graphID, graphIBCD, procDict)

    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    fasts._computePT_mut(zones, metrics, hook1)

    if tmy is None: return (t, tc, metrics)
    else: return (t, tc, metrics, tmy)

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
    if timelevel_target == timelevel_perfile or tc is None: 

       rank = Cmpi.rank
       tmp  = No_period*timelevel_period
       root = timelevel_perfile + ( (timelevel_motion - tmp)//timelevel_perfile)*timelevel_perfile

       FILE = tc_steady
       if os.access(FILE, os.F_OK): 
          tc = Cmpi.convertFile2SkeletonTree(FILE)
          tc = Cmpi.readZones(tc, FILE, rank=rank)

       FILE = directory+'/tc_'+str(root)+'.cgns'
       if os.access(FILE, os.F_OK): 
          tc_inst = Cmpi.convertFile2SkeletonTree(FILE)
          tc_inst = Cmpi.readZones(tc_inst, FILE, rank=rank)

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

       graphID   = Cmpi.computeGraph(tc, type='ID')
       graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
       procDict  = D2.getProcDict(tc)
       procList  = D2.getProcList(tc)
       graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }

       tc = Cmpi.convert2PartialTree(tc, rank=rank)

       #
       #Compactage tc
       # 
       # Get omp_mode
       omp_mode = 0
       node = Internal.getNodeFromName2(t, '.Solver#define')
       if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: omp_mode = Internal.getValue(node)

       # Reordone les zones pour garantir meme ordre entre t et tc
       FastI._reorder(t, tc, omp_mode)

       # Compactage arbre transfert
       zones = Internal.getZones(t)
       g     = graph['procDict']
       l     = graph['procList']
       
       X.miseAPlatDonorTree__(zones, tc, graph=graph)

       #Remise zero target
       if dtloc is not None: dtloc[4] = 0
       
    #
    #timelevel_motion larger than number of timelevels for 360degre 
    #
    if timelevel_motion > timelevel_360: dtloc[3] = 0  # remise a zero du compteur si 360degres 

    return tc, graph
