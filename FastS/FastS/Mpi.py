# FastS + MPI
from . import PyTree
from . import fasts
from .PyTree import display_temporal_criteria, createConvergenceHistory, extractConvergenceHistory, createStressNodes, createStatNodes, _computeStats, initStats, _computeEnstrophy, _computeVariables, _computeGrad, _compact, _applyBC, _init_metric, allocate_metric, _movegrid, _computeVelocityAle, copy_velocity_ale,  checkBalance, itt, distributeThreads, allocate_ssor, setIBCData_zero, display_cpu_efficiency
import timeit
import time as Time
import numpy

try:
    import Converter.PyTree as C
    import Converter.Mpi    as Cmpi
    import Distributor2.PyTree as D2
    import Converter.Internal as Internal
    import Connector.Mpi as Xmpi
    import Connector.PyTree as X
    import Connector
    import FastC.PyTree as FastC
    import os
    import math
except:
    raise ImportError("FastS: requires Converter, Connector and Distributor2 modules.")

try: range = xrange
except: pass

#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c", NIT=1, ucData=None):

    if isinstance(graph, list):
        #test pour savoir si graph est une liste de dictionnaires (explicite local)
        #ou juste un dictionnaire (explicite global, implicite)
        grapheliste=True
    else:
        grapheliste=False
    
    if graph is not None and grapheliste==False:
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
         distrib_omp = 0
         hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )
         nidom_loc = hook1["nidom_tot"]

         skip = 0
         if hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1: skip = 1

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
            #print('t_compute, rank= ', Time.time() - t0, Cmpi.rank)


            # dtloc GJeanmasson
            if exploc==1 and tc is not None:

               no_transfert = 1#comm_P2P
               process = Cmpi.rank

               zonesD    = Internal.getZones(tc)
               
               fasts.dtlocal2para_mpi(zones,zonesD,param_int_tc,param_real_tc,hook1,0,nstep,ompmode,no_transfert,process)

               if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
               elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 

               _applyBC(zones,metrics, hook1, nstep, ompmode, var=vars[0])    

               FastC.switchPointers2__(zones,nitmax,nstep)
                
               # Ghostcell

               if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
               elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1']
               timelevel_target = int(dtloc[4])

               if graph is not None and grapheliste == True:
                       procDict  = graph[nstep-1]['procDict']
                       graphID   = graph[nstep-1]['graphID']
                       graphIBCD = graph[nstep-1]['graphIBCD']
               else: 
                       procDict=None; graphID=None; graphIBCD=None          
                   
               _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict, nitmax, rk, exploc)

               no_transfert = 1#comm_P2P      

               fasts.recup3para_mpi(zones,zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1) 

               if nstep%2 == 0:
                   timelevel_target = int(dtloc[4])
                   vars = ['Density'  ]
                   if graph is not None and grapheliste == True:
                       procDict  = graph[nitmax+nstep-1]['procDict']
                       graphID   = graph[nitmax+nstep-1]['graphID']
                       graphIBCD = graph[nitmax+nstep-1]['graphIBCD']
                   else: 
                       procDict=None; graphID=None; graphIBCD=None      
                   _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, nitmax, hook1, graphID, graphIBCD, procDict, nitmax, rk, exploc, 2)

               if    nstep%2 == 0 and itypcp == 2: vars = ['Density'  ] 
               elif  nstep%2 == 1 and itypcp == 2: vars = ['Density_P1'] 
               _applyBC(zones,metrics, hook1, nstep, ompmode,  var=vars[0])



            else: ### Autres schemas
               # Ghostcell
               #tic=Time.time()  

               vars = FastC.varsP
               if nstep%2 == 0 and itypcp == 2: vars = FastC.varsN  # Choix du tableau pour application transfer et BC

               timelevel_target = int(dtloc[4])
              

               #print ('target= ', timelevel_target)
               #tic=Time.time()

               _fillGhostcells(zones, tc, metrics, timelevel_target , vars, nstep, ompmode, hook1, graphID, graphIBCD, procDict)

               # add unsteady Chimera transfers (motion) here
               if ucData is not None:
                   (graphX, intersectionDict, dictOfADT, 
                   dictOfNobOfRcvZones, dictOfNozOfRcvZones,
                   dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
                   dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC,
                   time, procDict, interpInDnrFrame, varType, tfreq) = ucData

                   if nstep%2 == 0 and itypcp == 2: 
                       VARS = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
                       if varType == 1: VARS += ['TurbulentSANuTilde']
                   else: 
                       VARS = ['Density_P1', 'VelocityX_P1', 'VelocityY_P1', 'VelocityZ_P1', 'Temperature_P1']
                       if varType == 1: VARS += ['TurbulentSANuTilde_P1']
                   for v in VARS: C._cpVars(t, 'centers:'+v, tc, v)
                   C._cpVars(t, "centers:cellN", tc, "cellN")

                   if nstep == 0 or nstep == nitmax or nstep%tfreq == 0:
                    Xmpi._transfer2(t, tc, VARS, graphX, intersectionDict, dictOfADT, 
                                    dictOfNobOfRcvZones, dictOfNozOfRcvZones,
                                    dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
                                    dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC, 
                                    time=time, absFrame=True,
                                    procDict=procDict, cellNName='cellN#Motion', 
                                    interpInDnrFrame=interpInDnrFrame, hook=hookTransfer)
               #print('t_transferts, rank= ', Time.time() - t0, Cmpi.rank)
               #toc1_=Time.time()-tic


      dtloc[3] +=1    #time_level_motion
      dtloc[4] +=1    #time_level_target

    else: 
      nstep_deb = 1
      nstep_fin = nitmax
      layer_mode= 1
      nit_c     = NIT
      FastC.HOOK["mpi"] = 1
      fasts._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, FastC.HOOK)

    #switch pointer a la fin du pas de temps
    if exploc==1 and tc is not None :
    #if exploc==2 and tc is not None and rk==3:
         if layer == 'Python':
             FastC.switchPointers__(zones, 1, 3)
         else:
             FastC.switchPointers3__(zones,nitmax)
    else:
         case = NIT%3
         #case = 2
         if case != 0 and itypcp < 2: FastC.switchPointers__(zones, case)
         if case != 0 and itypcp ==2: FastC.switchPointers__(zones, case, nitmax%2+2)

    # Difference avec sequentiel: pourquoi??
    #case = NIT%3
    #if case != 0: FastC.switchPointers__(zones, case, order=orderRk)

    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT          = 1
    return None
    #return tps_calcul, tps_com_transferts

#==============================================================================
def _fillGhostcells(zones, tc, metrics, timelevel_target, vars, nstep, omp_mode, hook1, graphID, graphIBCD, procDict, nitmax=1, rk=1, exploc=0, num_passage=1): 

   rank=Cmpi.rank

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

              for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              #if rank == 0: print('fillGC: timeleveltarget= ', timelevel_target)

              #recuperation Nb pas instationnaire dans tc
              type_transfert = 1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1,
                                        graph=graphIBCD, procDict=procDict)
              type_transfert = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, timelevel_target,#timecount,
                                        nstep, nitmax, rk, exploc, num_passage, varType=varType, compact=1, graph=graphID, procDict=procDict)

              #toc = Time.time() - tic
       # if rank == 0:
       #     print "Time in MPI send_buffer, irecv ","%.6f"%timecount[0]
       #     print "Time InterpTransfert (Inter)  ","%.6f"%timecount[1]
       #     print "Time InterpTransfert (Intra)  ","%.6f"%timecount[2]
       #     print "Time in getTransfersInter ","%.6f"%timecount[3]
       # if rank == 0: t0=timeit.default_timer()
       #apply BC
       #tic = Time.time()
       if exploc != 1:
           _applyBC(zones, metrics, hook1, nstep, omp_mode, var=vars[0])
       #toc = Time.time() - tic
       # if rank == 0:
       #     t1=timeit.default_timer()
       #     print "time/it (s) BC only (=t_RK X 1.0): ",(t1-t0)
   #return toc
   return None
#==============================================================================
def warmup(t, tc, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None):

    # Get omp_mode
    ompmode = PyTree.OMP_MODE
    node = Internal.getNodeFromName2(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastC._reorder(t, tc, ompmode)

    # Construction param_int et param_real des zones
    FastC._buildOwnData(t, Padding)

    
    # determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)

    
    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)# tab numpy
    nssiter = dtloc[0]   

    if isinstance(graph, list):
        #test pour savoir si graph est une liste de dictionnaires (explicite local)
        #ou juste un dictionnaire (explicite global, implicite)
        grapheliste=True
    else:
        grapheliste=False
    
    if graph is not None and grapheliste==False:
        procDict  = graph['procDict']
        graphID   = graph['graphID']
        graphIBCD = graph['graphIBCD']
    elif graph is not None and grapheliste==True:  ### Dans warmup tous les transferts doivent etre faits
        procDict  = graph[nssiter-2]['procDict']   ### On va chercher le graphe a nssiter-2 car a cette ssite tous les transferts
        graphID   = graph[nssiter-2]['graphID']    ### sont faits pour le schema a pas de temps local
        graphIBCD = graph[nssiter-2]['graphIBCD']
    else: 
        procDict=None; graphID=None; graphIBCD=None

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
    for nstep in range(1, int(dtloc[0])+1):
        hook1 = FastC.HOOK.copy()
        distrib_omp = 1
        hook1.update(  fasts.souszones_list(zones, metrics, FastC.HOOK, 1, nstep, distrib_omp) )

    _init_metric(t, metrics, ompmode)

    ssors = allocate_ssor(t, metrics, hook1, ompmode)

    # compact + align + init numa
    rmConsVars=True
    adjoint   =Adjoint

    t, FIRST_IT, zones2compact = FastC.createPrimVars(t, ompmode, rmConsVars, adjoint)
    FastC.HOOK['FIRST_IT']= FIRST_IT

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
        else: ale =True
        c += 1

    #
    # mise a jour vitesse entrainememnt
    #
    if ale == True and infos_ale is not None:
        print("ale actif. Teta et tetap=", infos_ale)
        teta = infos_ale[0];  tetap = infos_ale[1]
        FastC._motionlaw(t, teta, tetap)
        _computeVelocityAle(t,metrics)
    #
    # Compactage arbre transfert
    #
    if tc is not None:      
       X.miseAPlatDonorTree__(zones, tc, graph=graph,list_graph=list_graph)

       FastC.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int')[1]
       param_real_tc               = Internal.getNodeFromName1( tc, 'Parameter_real')
       if param_real_tc is not None: FastC.HOOK['param_real_tc']= param_real_tc[1]
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
        var = Internal.getNodesFromType1(sol[0] , 'DataArray_t')
        varmy=[]
        for v in var: varmy.append('centers:'+v[0])
        _compact(tmy, fields=varmy)

    #
    # remplissage ghostcells
    #
    hook1['lexit_lu'] = 0
    nstep             = 1
    nitrun            = 0
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4])

    _fillGhostcells(zones, tc, metrics, timelevel_target, ['Density'], nstep, ompmode, hook1,graphID, graphIBCD, procDict)


    #sys.exit()    

    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    fasts._computePT_mut(zones, metrics, hook1)

    return (t, tc, metrics)

#==============================================================================
# computeStress avec reduction 
#==============================================================================
def _computeStress(t, teff, metrics, xyz_ref=(0.,0.,0.)):
    """Compute efforts in teff.""" 
    ret = PyTree._computeStress(t, teff, metrics, xyz_ref)
    ret = numpy.array(ret, dtype=numpy.float64)
    ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    Cmpi.Allreduce(ret, ret1, Cmpi.SUM)
    return ret1.tolist()

#==============================================================================
# For periodic unsteady chimera join, parameter must be updated peridicaly 
#==============================================================================
def _UpdateUnsteadyJoinParam(t, tc, omega, timelevelInfos, graph, tc_steady='tc_steady.cgns', directory='.', iteration=0):

    bases = Internal.getNodesFromType1(t      , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    #print('own= ', own)
 
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


    if iteration == 0:    
       Node = Internal.getNodeFromName1(t, 'Iteration')
       if Node is not None: timelevel_motion = Internal.getNodeFromName1(t, 'Iteration')[1][0] - 1
       else: timelevel_motion = 0
    else:
       timelevel_motion = iteration

    timelevel_period = timelevelInfos["TimeLevelPeriod"]
    timelevel_360    = timelevelInfos["TimeLevel360"]
    timelevel_perfile= timelevelInfos["TimeLevelPerFile"]
    timelevel_axeRot = timelevelInfos["TimeLevelRotationAxis"]

    No_period = timelevel_motion//timelevel_period 
    #
    #target in no more in tc; need need data in a new file
    #

    if timelevel_target == timelevel_perfile or tc is None or timelevel_motion%timelevel_period == 0 : 

       #print('timelevel_target=',timelevel_target)
       rank = Cmpi.rank
       tmp  = No_period*timelevel_period
       root = timelevel_perfile + ( (timelevel_motion - tmp)//timelevel_perfile)*timelevel_perfile
       if root > timelevel_period : root=timelevel_period ### Comme 8000 pas multiple de 60 on force le load de tc_8000
       #print(root)
       if rank==0: print("timelevel_motion= ", timelevel_motion,"timelevel_target= ", timelevel_target)
    
       
       FILE = tc_steady
       if os.access(FILE, os.F_OK): 
          tc = Cmpi.convertFile2SkeletonTree(FILE)
          tc = Cmpi.readZones(tc, FILE, rank=rank)

       FILE = directory+'/tc_'+str(root)+'.cgns'
       if os.access(FILE, os.F_OK): 
          tc_inst = Cmpi.convertFile2SkeletonTree(FILE)
          tc_inst = Cmpi.readZones(tc_inst, FILE, rank=rank)

       if dtloc is not None:
           dtloc[4]=dtloc[4]%timelevel_perfile
       if timelevel_motion%timelevel_period == 0: ### On force dtloc[4] a 0 a chaque fin de periode azymutale
           dtloc[4]=0

       #
       #timelevel_motion larger than calculated peridicity; need to modify angle of rotation for azymuth periodicity
       #
       if timelevel_motion >= timelevel_period:
            No_period = timelevel_motion//timelevel_period
            iteration = timelevel_motion - No_period*timelevel_period 
       else:
            iteration = timelevel_motion


       if timelevel_motion >= timelevel_period: 
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
       
       tc = Internal.merge( [tc, tc_inst] )


       graphID   = Cmpi.computeGraph(tc, type='ID', it=iteration)
       graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
       procDict  = D2.getProcDict(tc)

       procList  = D2.getProcList(tc, sort=True)

       if rank==0 :  print ('proclist= ', procList)
   

       graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }


       tc = Cmpi.convert2PartialTree(tc, rank=rank)
 

       #Compactage tc
       # 
       # Get omp_mode
       ompmode = PyTree.OMP_MODE
       node = Internal.getNodeFromName2(t, '.Solver#define')
       if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

       # Reordone les zones pour garantir meme ordre entre t et tc
       FastC._reorder(t, tc, ompmode)

       # Compactage arbre transfert
       zones = Internal.getZones(t)

       #Remise zero target
       #if dtloc is not None: dtloc[4] = 0
       
       X.miseAPlatDonorTree__(zones, tc, graph=graph)

       #Remise zero target
       #if dtloc is not None: dtloc[4] = 0
       
       #timelevel_motion larger than number of timelevels for 360degre 
       #

    else:
       rank = Cmpi.rank


       if timelevel_motion >= timelevel_period:
            No_period = timelevel_motion//timelevel_period
            iteration = timelevel_motion - No_period*timelevel_period 
       else:
            iteration = timelevel_motion
    
       if rank==0:print('calcul du graph it', iteration)

         
       graphID   = Cmpi.computeGraph(tc, type='ID', it=iteration)
       graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
       procDict  = D2.getProcDict(tc)
       procList  = D2.getProcList(tc, sort=True)

       graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }




    if timelevel_motion > timelevel_360:
       print('remise a ZERO dans updateUnsteady')
       dtloc[3] = 0  # remise a zero du compteur si 360degres 

    return tc, graph

#====================================================================================
# Calcule la CFL en chaque maille et la met dans l'arbre t pour dtloc instationnaire
#====================================================================================
def computeCFL_dtlocal(t):

    import math
    import Transform.PyTree as T        
    from mpi4py import MPI  
 
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

               
    print(liste_BCPeriodiques)
                 

    (t,tc,metrics) = warmup(t,tc=None)

    zones = Internal.getZones(t)

    dim = Internal.getZoneDim(zones[0])
    
    print('dimPb= ', dimPb)
    
    fasts.prep_cfl(zones, metrics,1,1,1)

    #t = Internal.rmGhostCells(t, t, 2)
    """
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
    """

    cflmax_glob = 0.0
    cflmin_glob = 1.0e15 

    for z in zones:
        dim = Internal.getZoneDim(z)
        zp = T.subzone(z, (3,3,3), (dim[1]-2,dim[2]-2,dim[3]-2))
        cflmax = C.getMaxValue(zp, 'centers:CFL')
        cflmin = C.getMinValue(zp, 'centers:CFL')
        if cflmax > cflmax_glob:cflmax_glob = cflmax
        if cflmin < cflmin_glob:cflmin_glob = cflmin

    dtmin = 1.0/cflmax
    dtmax = 1.0/cflmin

    dts = numpy.array([dtmin,dtmax])
    dts_rcv = numpy.empty(6, dtype='d')

    MPI.COMM_WORLD.Allgather(dts,dts_rcv)

    dtmax = max(dts_rcv)
    dtmin = min(dts_rcv)
    

    timeLevel_max = math.floor(math.log((dtmax/dtmin))/math.log(2))

    print('dtmin, dtmax, Maxlevel = ', dtmin,',',dtmax,',',timeLevel_max)

    return (t,dtmin,timeLevel_max)

#==============================================================================
# decoupe maillage pour dtloc instationnaire
#==============================================================================
def _decoupe2(t, niveauMax, dtmin):

    import Transform.PyTree as T

    fes = Internal.getNodeFromName2(t,'FlowEquationSet')
    rs = Internal.getNodeFromName2(t,'ReferenceState')

    dtmin = 100000.
    taille_bloc = 25
    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t') 
    zones = []
    zones_decoupe = []

    dicoTps = {}
    numZone = []

    """
    liste_BCPeriodiques = []
    
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
         
    print(liste_BCPeriodiques)
    """     

    for b in bases:
        zones += Internal.getNodesFromType1(b, 'Zone_t')

    somme = 0
    numZone.append(somme)
    
    for z in zones:

        C._addVars(z,'centers:niveaux_temps')

        dim = Internal.getZoneDim(z)
        ni = dim[1]-5
        nj = dim[2]-5
        nk = dim[3]-5

        dimPb = dim[4]
        
        if dimPb==2:nk=1

        nbbloc_i = ni/taille_bloc
        nbbloc_i = int(nbbloc_i)
        if (nbbloc_i*taille_bloc < ni) : nbbloc_i = nbbloc_i + 1
        nbbloc_j = nj/taille_bloc
        nbbloc_j = int(nbbloc_j)
        if (nbbloc_j*taille_bloc < nj) : nbbloc_j = nbbloc_j + 1
        nbbloc_k = nk/taille_bloc
        nbbloc_k = int(nbbloc_k)
        if (nbbloc_k*taille_bloc < nk) : nbbloc_k = nbbloc_k + 1

        somme += nbbloc_i*nbbloc_j*nbbloc_k
        
        numZone.append(somme)

        print('decoupage de la zone', z[0])

        print(nbbloc_i,nbbloc_j,nbbloc_k)
        
         
        for k in range(0,nbbloc_k):
            for j in range(0,nbbloc_j):
                for i in range(0,nbbloc_i):

                    imin = i*taille_bloc + 3
                    jmin = j*taille_bloc + 3
                    kmin = k*taille_bloc + 3

                    imax = min((i+1)*taille_bloc+1,ni)
                    jmax = min((j+1)*taille_bloc+1,nj)
                    kmax = min((k+1)*taille_bloc+1,nk)

                    if dimPb==2:kmax=1
                   
                    zp = T.subzone(z,(imin,jmin,kmin),(imax,jmax,kmax))

                    cflmax_loc = C.getMaxValue(zp, 'centers:CFL')
                    dtmin_loc = 1./cflmax_loc
                    

                    niveauTps = math.log(((2**niveauMax)*dtmin)/dtmin_loc)/math.log(2)
     
 
                    fs = Internal.getNodeFromName(t,'FlowSolution#Centers')
                    nt = Internal.getNodeFromName(fs,'niveaux_temps')[1]
                    for i_ in range(imin, imax):
                        for j_ in range(jmin, jmax):
                            for k_ in range(kmin, kmax):
                                 nt[i,j,k] = float(niveauTps)


    return t
