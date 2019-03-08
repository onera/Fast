# FastS + MPI
import PyTree
import fasts
from PyTree import display_temporal_criteria, createConvergenceHistory, extractConvergenceHistory, createStressNodes, _computeStress, createPrimVars, _createPrimVars, createStatNodes, _computeStats, initStats, _computeEnstrophy, _computeVariables, _computeGrad, _compact, _applyBC, _buildOwnData, _init_metric, allocate_metric, _BCcompact, _motionlaw, _movegrid, _computeVelocityAle, checkBalance, itt, HOOK, distributeThreads, _build_omp, allocate_ssor, setIBCData_zero, display_cpu_efficiency
import timeit
import time as Time
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
    raise ImportError("FastS: requires Converter, Connector and Distributor2 modules.")

#==============================================================================
#def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c", NIT=1, tps_calcul=None, tps_com_transferts=None):
def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c", NIT=1):
    if graph is not None:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
    else: 
        procDict=None; graphID=None; graphIBCD=None

    base = Internal.getNodeFromType1(t   ,"CGNSBase_t")
    own  = Internal.getNodeFromName1(base, '.Solver#ownData')  
    dtloc= Internal.getNodeFromName1(own , '.Solver#dtloc')

    zones= Internal.getZones(t)

    node = Internal.getNodeFromName1(base, '.Solver#define')
    omp_node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode  = PyTree.OMP_MODE
    if  omp_node is not None: ompmode = Internal.getValue(omp_node)

    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])


    #### a blinder...
    param_int_firstZone = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1]
    itypcp = param_int_firstZone[29]
    rk     = param_int_firstZone[52]
    exploc = param_int_firstZone[54]
    #### a blinder...

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

      for nstep in xrange(1, nitmax+1): # pas RK ou ssiterations


         hook1 = PyTree.HOOK.copy()
         distrib_omp = 0
         hook1.update(  fasts.souszones_list(zones, metrics, PyTree.HOOK, nitrun, nstep, distrib_omp) )
         nidom_loc = hook1["nidom_tot"]

         skip = 0
         if (hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1): skip = 1

         # calcul Navier Stokes + appli CL
         if nidom_loc > 0 and skip ==0:
            # Navier-Stokes
            nstep_deb = nstep
            nstep_fin = nstep
            layer_mode= 0
            nit_c     = 1
            #t0=Time.time()
            #tic=Time.time()
            fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, hook1)
            #toc1=Time.time()-tic  
            #print 't_compute, rank= ', Time.time() - t0, Cmpi.rank


            #dtloc GJeanmasson
            if exploc==2 and tc is not None and rk==3:
               fasts.dtlocal2para_(zones, zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1, dest)

               if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
               elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
               _applyBC(zones,metrics, hook1, nstep, ompmode, var=vars[0], rk=rk, exploc=exploc)    

               FastI.switchPointers2__(zones,nitmax,nstep)
                
               # Ghostcell
               if (nitmax%3 != 0) : # Tous les schemas sauf constantinescu RK3
                   if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
                   elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1']
                   timelevel_target = int(dtloc[4])
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict, nitmax, rk, exploc)
            
               fasts.recup3para_(zones,zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1) 

               if (nstep%2==0) :
                   timelevel_target = int(dtloc[4])
                   vars = ['Density'  ]
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, nitmax, hook1, graphID, graphIBCD, procDict, nitmax, rk, exploc, 2)

               if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
               elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
               _applyBC_(zones,metrics, hook1, nstep, ompmode,  var=vars[0])



            else:
               # Ghostcell
               #tic=Time.time()  
               vars = PyTree.varsP
               if  nstep == 2 and itypcp == 2 : vars = PyTree.varsN  # Choix du tableau pour application transfer et BC
               timelevel_target = int(dtloc[4])
               tic=Time.time()
               _fillGhostcells(zones, tc, metrics, timelevel_target , vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict)
               #print 't_transferts, rank= ', Time.time() - t0, Cmpi.rank
               #toc1_=Time.time()-tic  
                        
               #tps_calcul = tps_calcul + toc1 + toc3  
               #tps_com_transferts = tps_com_transferts + toc1_ - toc3 

    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT
      PyTree.HOOK["mpi"] = 1
      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, PyTree.HOOK)

    #switch pointer a la fin du pas de temps
    if exploc==2 and tc is not None and rk==3:
         if layer == 'Python':
             FastI.switchPointers__(zones, 1, 3)
         else:
             FastI.switchPointers3__(zones,nitmax)
    else:
         case = NIT%3
         #case = 2
         if case != 0 and itypcp < 2: FastI.switchPointers__(zones, case)
         if case != 0 and itypcp ==2: FastI.switchPointers__(zones, case, nitmax%2+2)

    # Difference avec sequentiel: pourquoi??
    #case = NIT%3
    #if case != 0: FastI.switchPointers__(zones, case, order=orderRk)

    # flag pour derivee temporelle 1er pas de temps implicit
    PyTree.HOOK["FIRST_IT"]  = 1
    PyTree.FIRST_IT          = 1
    return None
    #return tps_calcul, tps_com_transferts


#==============================================================================
#def _compute(t, metrics, nitrun, tc=None, graph=None, layer="Python", NIT=1):
def _computeguillaume1(t, metrics, nitrun, tc=None, graph=None, layer="Python", NIT=1, list_graph=None, tps_calcul=None, tps_com_transferts=None):
    #if graph is not None:
    #    procDict  = graph['procDict']
    #    graphID   = graph['graphID']
    #    graphIBCD = graph['graphIBCD']
    #else: 
    #    procDict=None; graphID=None; graphIBCD=None

    base = Internal.getNodeFromType1(t,"CGNSBase_t")
    own  = Internal.getNodeFromName1(base, '.Solver#ownData') 
    dtloc= Internal.getNodeFromName1(own, '.Solver#dtloc')
    zones= Internal.getZones(t)
    node = Internal.getNodeFromName(t, '.Solver#define')
    node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode = PyTree.OMP_MODE
    if node is not None: ompmode = Internal.getValue(node)
    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])
    orderRk = int(dtloc[len(dtloc)-1])

    node = Internal.getNodeFromName1(base, '.Solver#define')
    exp_local = Internal.getNodeFromName1(node, 'exp_local') 
    explocal = Internal.getValue(exp_local)
    r_k = Internal.getNodeFromName1(node, 'rk') 
    rk = Internal.getValue(r_k)

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

    if (layer=='Python'):
        taille_tabs = 2000000000
        rostk = numpy.empty(taille_tabs, dtype=numpy.float64)  # tableau de stockage des valeurs 
        drodmstk = numpy.empty(taille_tabs, dtype=numpy.float64) # tableau de stockage des flux (drodm)
        constk = numpy.empty(taille_tabs, dtype=numpy.float64) # tableau de stockage des flux (conservativite

    zonesD    = Internal.getZones(tc)

    nitmaxs = 1

    if layer == "Python": 

 
      #tps_calcul = 0.0
      #tps_com_transferts = 0.0 
        
      for nstep in xrange(1, nitmax+1): # pas RK ou ssiterations

         #print 'coucou', Cmpi.rank

         if graph is not None:
             procDict  = list_graph[nstep-1]['procDict']
             graphID   = list_graph[nstep-1]['graphID']
             graphIBCD = list_graph[nstep-1]['graphIBCD']
         else: 
             procDict=None; graphID=None; graphIBCD=None          

         #print 'nstep= ', nstep

         hook1 = PyTree.HOOK.copy()
         distrib_omp = 0
         hook1.update(  fasts.souszones_list(zones, metrics, PyTree.HOOK, nitrun, nstep, distrib_omp) )
         nidom_loc = hook1["nidom_tot"]

         skip = 0
         if (hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1): skip = 1

         # calcul Navier Stokes + appli CL
         if nidom_loc > 0 and skip ==0:
            # Navier-Stokes
            nstep_deb = nstep
            nstep_fin = nstep
            layer_mode= 0
            nit_c     = 1
            
            #tic=timeit.default_timer()
            #tic=Time.time()
            fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, hook1)
            #toc1=timeit.default_timer()-tic 
            #toc1=Time.time()-tic           
 

            #tic=Time.time()
            if    (nstep%2 == 0)  and itypcp == 2 : vars = PyTree.varsN 
            elif  (nstep%2 == 1)  and itypcp == 2 : vars = PyTree.varsP 
            _applyBC(zones,metrics, hook1, nstep, ompmode, var=vars[0], rk=rk, exploc=explocal)
            #toc2=Time.time()-tic              
            #print 'coucou'

            #tic=Time.time()
            #for comm_P2P in xrange(1,param_int[0]+1):
            no_transfert = 1#comm_P2P
            process = Cmpi.rank
            #print Cmpi.rank, comm_P2P
            #tic=Time.time()
            fasts.dtlocal2para_mpi(zones,zonesD,param_int,param_real,hook1,rostk,drodmstk,constk,0,nstep,ompmode,taille_tabs,no_transfert,process)
          
            #print 'coucou'
            FastI.switchPointers2__(zones,nitmax,nstep)
            #toc3=Time.time()-tic
                
            # Ghostcell
            #tic=Time.time()
            if (nstep != nitmax) : # Tous les schemas sauf constantinescu RK3
                 if    (nstep%2 == 0)  and itypcp == 2 : vars = PyTree.varsN  # Choix du tableau pour application transfer et BC
                 elif  (nstep%2 == 1)  and itypcp == 2 : vars = PyTree.varsP
                 timelevel_target = int(dtloc[4])
                 _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict, nitmax, rk, explocal)
            #toc1_=Time.time()-tic
            

            #tic=Time.time()
            #for comm_P2P in xrange(1,param_int[0]+1):
            no_transfert = 1#comm_P2P           
            fasts.recup3para_mpi(zones,zonesD,param_int,param_real,hook1,rostk,0,nstep,ompmode,taille_tabs,no_transfert) 
            #toc4=Time.time()-tic

            #tic=Time.time()
            if (nstep%2==0) :
                timelevel_target = int(dtloc[4])
                vars = PyTree.varsN
                if graph is not None:
                    procDict  = list_graph[nitmax+nstep-1]['procDict']
                    graphID   = list_graph[nitmax+nstep-1]['graphID']
                    graphIBCD = list_graph[nitmax+nstep-1]['graphIBCD']
                else: 
                    procDict=None; graphID=None; graphIBCD=None      
                _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, nitmax, hook1, graphID, graphIBCD, procDict, nitmax, rk, explocal, 2)
            #toc2_=Time.time()-tic 
            #print 'coucou'

            #tps_calcul = tps_calcul + toc1 + toc2 + toc4  + toc3 
            #tps_com_transferts = tps_com_transferts + toc1 #toc2_+ toc1_ 


    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT
      PyTree.HOOK["mpi"] = 1
      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, PyTree.HOOK)

      FastI.switchPointers3__(zones,nitmax)
      #FastI.switchPointers__(zones, 1, 3)                       
   
    if (layer == 'Python') :
      #FastI.switchPointers3__(zones,nitmax)
      FastI.switchPointers__(zones, 1, 3)

    # switch pointers
    #case = NIT%3
    #tic=Time.time()
    #if case != 0: FastI.switchPointers__(zones,1,3)
    #toc5=Time.time()-tic
    #tps_calcul = tps_calcul + toc5
    # flag pour derivee temporelle 1er pas de temps implicit
    PyTree.HOOK["FIRST_IT"]  = 1
    PyTree.FIRST_IT          = 1

    #print 'rang, tps calcul : ', Cmpi.rank, tps_calcul
    #print 'rang, tps transfert : ', Cmpi.rank, tps_com_transferts
    return tps_calcul, tps_com_transferts
    #return None
#==============================================================================


def _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, omp_mode, hook1, graphID, graphIBCD, procDict, nitmax=1, rk=1, exploc=1, num_passage=1): 

   rank=Cmpi.rank 

   #print vars

   #timecount = numpy.zeros(4, dtype=numpy.float64)
   timecount = []

   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None:
           tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
           if tc_compact is not None:
              param_real= tc_compact[1]              
              param_int = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1['neq_max'] == 5: varType = 2
              else                    : varType = 21

              bcType = PyTree.HOOKIBC[0]; Gamma=PyTree.HOOKIBC[1]; Cv=PyTree.HOOKIBC[2]; Mus=PyTree.HOOKIBC[3]; Cs=PyTree.HOOKIBC[4]; Ts=PyTree.HOOKIBC[5]
       
              #tic = Time.time()              
              if (rk==3 and exploc==2):
                 if nstep <= nitmax: 
                      for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)
              else:
                 if nstep <= 3: 
                      for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)
              #toc = Time.time() - tic


              # #recuperation Nb pas instationnaire dans tc
              type_transfert = 1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        nstep, nitmax, rk, exploc, num_passage, bcType=bcType, varType=varType, compact=1,
                                        graph=graphIBCD, procDict=procDict, Gamma=Gamma, Cv=Cv, MuS=Mus, Cs=Cs, Ts=Ts)
              type_transfert = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphID, procDict=procDict)

              #toc = Time.time() - tic
       # if (rank == 0):
       #     print "Time in MPI send_buffer, irecv ","%.6f"%timecount[0]
       #     print "Time InterpTransfert (Inter)  ","%.6f"%timecount[1]
       #     print "Time InterpTransfert (Intra)  ","%.6f"%timecount[2]
       #     print "Time in getTransfersInter ","%.6f"%timecount[3]
       # if (rank == 0 ): t0=timeit.default_timer()
       #apply BC
       #tic = Time.time()
       if (exploc != 2):
           _applyBC(zones, metrics, hook1, nstep, omp_mode, var=vars[0], rk=1, exploc=1)
       #toc = Time.time() - tic
       # if (rank == 0 ):
       #     t1=timeit.default_timer()
       #     print "time/it (s) BC only (=t_RK X 1.0): ",(t1-t0)
   #return toc
   return None
#==============================================================================
def warmup(t, tc, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None):

    if graph is not None:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
    else: 
        procDict=None; graphID=None; graphIBCD=None

    # Get omp_mode
    ompmode = PyTree.OMP_MODE
    node = Internal.getNodeFromName2(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastI._reorder(t, tc, ompmode)

    # Construction param_int et param_real des zones
    _buildOwnData(t, Padding)

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
        distrib_omp = 1
        hook1.update(  fasts.souszones_list(zones, metrics, PyTree.HOOK, 1, nstep, distrib_omp) )

    _init_metric(t, metrics, hook1, ompmode)

    ssors = allocate_ssor(t, metrics, hook1, ompmode)

    # compact + align + init numa
    rmConsVars=True
    adjoint   =Adjoint
    t = createPrimVars(t, ompmode, rmConsVars, adjoint)

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
       X.miseAPlatDonorTree__(zones, tc, graph=graph,list_graph=list_graph)

       PyTree.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int')[1]
       param_real_tc               = Internal.getNodeFromName1( tc, 'Parameter_real')
       if param_real_tc is not None: PyTree.HOOK['param_real_tc']= param_real_tc[1]
    else:
        PyTree.HOOK['param_real_tc'] = None
        PyTree.HOOK['param_int_tc']  = None 

    if ssors is not []:
        PyTree.HOOK['ssors'] = ssors
    else:
        PyTree.HOOK['ssors'] = None
    #
    # remplissage ghostcells
    #
    hook1['lexit_lu'] = 0
    nstep             = 1
    nitrun            = 0
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4]) 

    _fillGhostcells(zones, tc, metrics, timelevel_target, ['Density'], nstep, ompmode, hook1,graphID, graphIBCD, procDict)

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
       ompmode = PyTree.OMP_MODE
       node = Internal.getNodeFromName2(t, '.Solver#define')
       if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

       # Reordone les zones pour garantir meme ordre entre t et tc
       FastI._reorder(t, tc, ompmode)

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

