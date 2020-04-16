"""Common functions for FAST solvers.
"""
from . import fast
from . import Fast
__version__ = Fast.__version__

OMP_MODE = 0

try: 
    import FastC.PyTree as FastC
    from FastC.PyTree import _setNum2Zones, _setNum2Base, load, save, loadFile, saveFile, loadTree, saveTree
     
except: 
    raise ImportError("Fast.PyTree: requires FastC module.")

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import FastP.PyTree as FastP
    import FastS.PyTree as FastS  
    import FastLBM.PyTree as FastLBM
    import Connector.PyTree as X
except:
    raise ImportError("Fast.PyTree: requires Converter and FastSP modules.")

import numpy

#==============================================================================
# Initialisation parametre calcul: calcul metric + var primitive + compactage 
# + alignement + placement DRAM
#==============================================================================
def warmup(t, tc=None, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None):
    """Perform necessary operations for the solver to run."""

    #global FastC.FIRST_IT, FastC.HOOK
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
    nssiter = int(dtloc[0])

    zones = Internal.getZones(t)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it

    # allocation d espace dans param_int pour stockage info openmp
    FastC._build_omp(t) 

    
    zones_unstr =[]; zones_str = []; zones_ns = []; zones_lbm = []
    for z in zones:
       ztype = Internal.getValue(  Internal.getNodeFromName(z, 'ZoneType') )
       if ztype=='Structured':
         zones_str.append(z)
         param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
         if param_int[27] == 4:   #IFLOW=4
            zones_lbm.append(z)
         else:
            zones_ns.append(z)
            
       else:
         zones_unstr.append(z)

    # alloue metric: tijk, ventijk, ssiter_loc
    # init         : ssiter_loc
    metrics = allocate_metric(t)

    metrics_str=[]; metrics_unstr=[]; metrics_ns = []; metrics_lbm = []
    c  =0
    for z in zones:
       ztype = Internal.getValue(  Internal.getNodeFromName(z, 'ZoneType') )
       if ztype=='Structured':
         metrics_str.append( metrics[c])
         param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
         if param_int[27] == 4:   #IFLOW=4
            metrics_lbm.append(z)
         else:
            metrics_ns.append(z)
       else:
         metrics_unstr.append( metrics[c])
       c+=1
    # Contruction BC_int et BC_real pour CLa
    FastC._BCcompact(zones_str) 
    FastC._BCcompactNG(zones_unstr) 

    #determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)
    #print 'int(dtloc[0])= ', int(dtloc[0])
    if len(zones_str) != 0: 
       for nstep in range(1, int(dtloc[0])+1):
           hook1       = FastC.HOOK.copy()
           distrib_omp = 1
           hook1.update(  FastS.fasts.souszones_list(zones_str, metrics_str, FastC.HOOK, 1, nstep, distrib_omp) )   
    
    #init metric 
    FastS.fasts.init_metric(zones_str  , metrics_str  , ompmode)
    FastP.fastp.init_metric(zones_unstr, metrics_unstr)

    # compact + align + init numa
    rmConsVars=True
    adjoint=Adjoint

    t, FastC.FIRST_IT, zones2compact = FastC.createPrimVars(t, ompmode, rmConsVars, adjoint)
    FastC.HOOK['FIRST_IT']= FastC.FIRST_IT

    if len(zones_str) == 0: hook1 = FastC.HOOK.copy(); hook1.update( {'nidom_tot':len(zones), 'lexit_lu':0, 'lssiter_verif':0})

    #compactage des champs en fonction option de calcul  
    count = -1
    if ompmode == 1: count = 0          
    for data in zones2compact:
        if ompmode == 1: count += 1
        zone    = data[0]
        varnames= data[1]
        for fields in varnames:
            _compact(zone, fields=fields, mode=count)

    #on recupere les zones a nouveau car create primvar rend zones caduc
    zones = Internal.getZones(t)
    #zones_str, zones_unstr, metrics_str, metrics_unstr = tri_zones( zones, metrics)
    infos_zones = tri_zones( zones, metrics)  #info-zone: dico 4 item contenant la list [ zones, metrics] pour Struct, unstruc, LBM, NS
    # init Q variables from macro variable if lbm
    #FastLBM._init_Q(infos_zones["LBM"][0])
    FastLBM._init_Q(infos_zones["LBM"][0])

    #allocate tab ssor en structure
    zones_str   = infos_zones["struct"][0]
    metrics_str = infos_zones["struct"][1]
    ssors = FastS.fasts.allocate_ssor(zones_str, metrics_str, nssiter, hook1, ompmode)

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
        else: ale = True
        c += 1

    #
    # mise a jour vitesse entrainememnt
    #
    if ale == True and infos_ale is not None:
        print("ale actif. Teta et tetap=", infos_ale)
        teta = infos_ale[0];  tetap = infos_ale[1]
        FastC._motionlaw(t, teta, tetap)
        if len(zones_str)   != 0: FastS.fasts.computePT_velocity_ale(zones_str,  metrics_str, FastC.HOOK, ompmode)
        if len(zones_unstr) != 0: print("coder fastp.computePT_velocity_ale(zones_unstr,  metrics_unstr")
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
       FastS._createTBLESA(tc,nbpts_linelets)
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

    if   len( infos_zones["LBM"][0]) ==0:  vars=FastC.varsN + FastC.varsMLBM
    elif len( infos_zones["NS"][0]) ==0 and len( infos_zones["unstruct"][0]) ==0: vars=FastC.varsMLBM + FastC.varsN
    else:
      print("echange NS-LBM a coder")
      import sys; sys.exit()
    
    _fillGhostcells(zones, tc,  infos_zones, timelevel_target, vars, nstep,  ompmode, hook1) 
    if tc is not None: C._rmVars(tc, 'FlowSolution')

    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    FastS.fasts._computePT_mut(infos_zones['struct'][0], infos_zones['struct'][1], hook1)
    #FastP.fastp._computePT_mut(infos_zones['unstruct'][0], infos_zones['unstruct'][1], hook1)
    print("decommenter _computePT_mut dans warmup Fast pour zone non structure.")

    return (t, tc, metrics)

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

        ztype = Internal.getValue(  Internal.getNodeFromName(z, 'ZoneType') )
        if ztype=='Structured':
           param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]

           if param_int[27] == 4:   #IFLOW=4
              metrics.append( FastLBM.fastlbm.allocate_metric(z, nssiter))
           else:
              metrics.append( FastS.fasts.allocate_metric(z, nssiter))
        else:
           metrics.append(FastP.fastp.allocate_metric(z, nssiter))
    return metrics

#==============================================================================
def _compact(t, containers=[Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__], fields=None, mode=None, init=True):
    if  mode is not None:
      if mode == -1: thread_numa = 0
      else: thread_numa = 1
    else: thread_numa = 0

    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getValue(  Internal.getNodeFromName(z, 'ZoneType') )
        ars   = FastC.getFields2Compact__(z, containers, fields)
        sh    = None ; size = None
        val   = [] # valid fields
        for a in ars:
            a1 = a[1]
            if sh is None: sh = a1.shape; size = a1.size; val.append(a)
            elif a1.shape == sh: val.append(a)
        nfields = len(val)
        if nfields > 0:
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')  # noeud
            # Create an equivalent contiguous numpy [flat]
            eq = numpy.empty(nfields*(size+ param_int[1][66]), dtype=numpy.float64) # add a shift  between prim. variables (param_int[1][SHIFTVAR])
            c = 0        
            if param_int is None:
                raise ValueError("_compact: Parameter_int is missing for zone %s."%z[0])
            for a in val:
                a1 = a[1]
                # Copy elements
                ptr = a1.reshape((size), order='F') # no copy I hope

                if init: 
                  if ztype=='Structured':
                    FastS.fasts.initNuma( ptr, eq, param_int, c, thread_numa )
                  else:
                    #fastp.initNuma( ptr, eq, param_int, c, thread_numa )
                    eq[c*size:(c+1)*size] = ptr[:]

                # Replace numpys with slice
                a[1] = eq[c*(size)+c*param_int[1][66]:(c+1)*(size)+c*param_int[1][66]]
                a[1] = a[1].reshape(sh, order='F')

                c += 1
    return None

#==============================================================================
def _fillGhostcells(zones, tc, infos_zones, timelevel_target, vars, nstep, ompmode, hook1, nitmax=1, rk=1, exploc=1, num_passage=1): 

   if hook1['lexit_lu'] ==0:

       #transfert
       if tc is not None :
           tc_compact = Internal.getNodeFromName1( tc, 'Parameter_real')
           #Si param_real n'existe pas, alors pas de raccord dans tc
           if tc_compact is not  None:

              param_real= tc_compact[1]
              param_int = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1["neq_max"] == 5: varType =  2
              else                    : varType = 21
              if vars[0] == 'Q1_M1'   : varType =  4
              #for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)
              C._cpVars(zones, 'centers:'+vars[0], zonesD, vars[0])

              type_transfert = 2  # 0= ID uniquement, 1= IBC uniquement, 2= All
              no_transfert   = 1  # dans la list des transfert point a point
              X.connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage)#,timecount)
       #apply BC
       if rk != 3 and exploc != 2:
          _applyBC(infos_zones, hook1, nstep, ompmode, var= vars)    

            

   return None
#==============================================================================
# compute in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c", NIT=1):
    """Compute a given number of iterations."""

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own     , '.Solver#dtloc')    # noeud

    node = Internal.getNodeFromName1(bases[0], '.Solver#define')
    omp_node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode  = OMP_MODE
    if  omp_node is not None: ompmode = Internal.getValue(omp_node)

    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])                 

    zones       = Internal.getZones(t)
    infos_zones = tri_zones( zones, metrics)

    #### a blinder...
    param_int_firstZone = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1]
    iflow  = param_int_firstZone[27]
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

         hook1 = FastC.HOOK.copy()
         distrib_omp = 0
         #hook1.update(  FastS.fasts.souszones_list(zones_str, metrics_str, FastC.HOOK, nitrun, nstep, distrib_omp) )
         hook1.update(  FastS.fasts.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, distrib_omp) )

         #nidom_loc = hook1["nidom_tot"] + len(zones_unstr)
         nidom_loc = hook1["nidom_tot"]
         #print("nidom_loc=", nidom_loc)
        
         skip      = 0
         if hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1: skip = 1

         # calcul Navier_stokes + appli CL
         if nidom_loc > 0 and skip == 0:
            nstep_deb = nstep
            nstep_fin = nstep
            layer_mode= 0
            nit_c     = 1
            #t0=Time.time()
            fast._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c, hook1)
            #print('t_compute = %f'%(Time.time() - t0))

            # dtloc GJeanmasson
            if exploc==2 and tc is not None and rk==3:
               FastS.fasts.dtlocal2para_(zones, zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1, dest)

               if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ] 
               elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1'] 
               _applyBC(infos_zones, hook1, nstep, ompmode, var=vars[0])    

               FastC.switchPointers2__(zones_str,nitmax,nstep)
                
               # Ghostcell
               if nitmax%3 != 0: # Tous les schemas sauf constantinescu RK3
                   if    (nstep%2 == 0)  and itypcp == 2 : vars = ['Density'  ]  # Choix du tableau pour application transfer et BC
                   elif  (nstep%2 == 1)  and itypcp == 2 : vars = ['Density_P1']
                   timelevel_target = int(dtloc[4])
                   _fillGhostcells(zones, tc, infos_zones, timelevel_target, vars, nstep, ompmode, hook1, nitmax, rk, exploc)
            
               FastS.fasts.recup3para_(zones,zones_tc, param_int_tc, param_real_tc, hook1, 0, nstep, ompmode, 1) 

               if nstep%2 == 0:
                   timelevel_target = int(dtloc[4])
                   vars = ['Density'  ]
                   _fillGhostcells(zones, tc,  infos_zones, timelevel_target, vars, nstep, nitmax, hook1, nitmax, rk, exploc, 2)

               if    nstep%2 == 0 and itypcp == 2 : vars = ['Density'  ] 
               elif  nstep%2 == 1 and itypcp == 2 : vars = ['Density_P1'] 
               _applyBC(infos_zones, hook1, nstep, ompmode,  var=vars[0])


            else:
              #Ghostcell
              varsNS = FastC.varsP 
              if nstep%2 == 0 and itypcp == 2 : varsNS = FastC.varsN  # Choix du tableau pour application transfer et BC
              vars = varsNS + FastC.varsMLBM
              if iflow == 4: vars= FastC.varsMLBM + varsNS
              timelevel_target = int(dtloc[4])
              _fillGhostcells(zones, tc,  infos_zones, timelevel_target, vars, nstep, ompmode, hook1)

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

      fast._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, ompmode, nit_c , FastC.HOOK)

    zones_ns    = infos_zones["NS"][0]
    metrics_ns  = infos_zones["NS"][1]
    zones_lbm   = infos_zones["LBM"][0]
    metrics_lbm = infos_zones["LBM"][1]
    #switch pointer a la fin du pas de temps
    if exploc==2 and tc is not None and rk==3:
         if layer == 'Python':
             FastC.switchPointers__(zones_ns, 1, 3)
         else:
             FastC.switchPointers3__(zones_ns,nitmax)
    else:
         case = NIT%3
         if case != 0 and itypcp < 2: FastC.switchPointers__(zones_ns, case)
         if case != 0 and itypcp ==2: FastC.switchPointers__(zones_ns, case, nitmax%2+2)

    if  NIT%2 != 0:
      FastC.switchPointersLBM__(zones_lbm, FastLBM.NQ)

    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT          = 1

    return None

#==============================================================================
# applyBC
#==============================================================================
def _applyBC(infos_zones, hook1, nstep, ompmode, var=["Density","Q1"]):

    if 'Q' in var[0]: 
         varlbm =  var[0]
         varns  =  var[1]
    else: 
         varlbm =  var[1]
         varns  =  var[0]

    FastS.fasts._applyBC(infos_zones["struct"][0]  , infos_zones["struct"][1]    , hook1, nstep, ompmode, varns  )
    FastLBM.fastlbm._applyBC(infos_zones["LBM"][0] , infos_zones["LBM"][1]       , hook1, nstep, ompmode, varlbm )
    FastP.fastp._applyBC(infos_zones["unstruct"][0], infos_zones["unstruct"][1]  , hook1, nstep, ompmode, varns  )
    

    return None

#==============================================================================
# renvoi list zone structure et polyhedrique a partir arbre hybride
#==============================================================================
def tri_zones(zones, metrics):

    infos_zones = {}
    zones_unstr =[]; zones_str = []; zones_ns = []; zones_lbm = []
    metrics_unstr =[]; metrics_str = []; metrics_ns = []; metrics_lbm = []

    c=0
    for z in zones:
       ztype = Internal.getValue(  Internal.getNodeFromName(z, 'ZoneType') )
       if ztype=='Structured':
         zones_str.append(z)
         metrics_str.append( metrics[c])
         param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
         if param_int[27] == 4:   #IFLOW=4
            zones_lbm.append(z)
            metrics_lbm.append( metrics[c] )
         else:
            zones_ns.append(z)
            metrics_ns.append( metrics[c] )
       else:
         zones_unstr.append(z)
         metrics_unstr.append( metrics[c])
         zones_ns.append(z)
         metrics_ns.append( metrics[c] )
       c+=1
    infos_zones['struct'  ]=[zones_str, metrics_str]
    infos_zones['unstruct']=[zones_unstr, metrics_unstr]
    infos_zones['LBM'     ]=[zones_lbm, metrics_lbm]
    infos_zones['NS'      ]=[zones_ns , metrics_ns]
    return infos_zones
