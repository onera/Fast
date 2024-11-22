"""Common functions for FAST solvers.
"""
import os
import numpy

import FastC.fastc
import FastS.fasts

from . import fast
from . import Fast
__version__ = Fast.__version__

OMP_MODE = 0
try:
    import FastC.PyTree as FastC
    from FastC.PyTree import _setNum2Zones, _setNum2Base, setNum2Zones, setNum2Base, load, \
        save, loadFile, loadFileG, saveFile, loadTree, saveTree, \
        getDictOfNobNozOfRcvZones, _addPair, getDictOfNobNozOfDnrZones, _pushCenters, \
        ramp, _rampTimeStep, _rampCFL, _pushWalls
except ImportError:
    raise ImportError("Fast.PyTree: requires FastC module.")

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import FastS.PyTree as FastS
    #import FastP.PyTree as FastP
    import Connector.PyTree as X
    from . import VariablesSharePyTree as VSHARE
except:
    raise ImportError("Fast.PyTree: requires Converter, Connector, FastS, FastASLBM modules.")

FASTLBM = False
try:
    import FastASLBM.PyTree as FastLBM
    FASTLBM = True
except:
    FASTLBM = False

try:
    OMP_NUM_THREADS = os.environ['OMP_NUM_THREADS']
    OMP_NUM_THREADS = int(OMP_NUM_THREADS)
except: OMP_NUM_THREADS = 1

#==============================================================================
# Initialisation parametre calcul: calcul metric + var primitive + compactage
# + alignement + placement DRAM
#==============================================================================
def warmup(t, tc=None, graph=None, infos_ale=None, Adjoint=False, tmy=None, list_graph=None, Padding=None, SizeBlockTarget=1000, nghost=2, verbose=0):
    """Perform necessary operations for the solver to run."""

    # Barriere si version sans FastP
    for z in Internal.getZones(t):
        if Internal.getZoneType(z) == 2:
            raise TypeError('warmup: unstructured zones detected and no FastP. Please suppress them from pyTree.')

    # renumerotation cellule/face si necessaire pour HPC (threading et cache blocking)
    node = Internal.getNodeFromName(t, 'HpcSplitInfo')
    lsplit = True
    if node is not None:
        Nthread = node[1][0]
        Ncache  = node[1][1]
        if Ncache == SizeBlockTarget and Nthread == OMP_NUM_THREADS: lsplit = False
        print("INFOSPLIT initial: Nthread=",Nthread,'Ncache=',Ncache,'NewSplit=',lsplit)

    #if lsplit: t,tc = FastP.SplitThreadCache(t, tc, SizeBlockTarget)

    #global FastC.FIRST_IT, FastC.HOOK
    # Get omp_mode
    ompmode = OMP_MODE
    node = Internal.getNodeFromName1(t, '.Solver#define')
    if node is not None:
        node = Internal.getNodeFromName1(node, 'omp_mode')
        if  node is not None: ompmode = Internal.getValue(node)

    # Reordone les zones pour garantir meme ordre entre t et tc
    FastC._reorder(t, tc)

    # Construction param_int et param_real des zones
    FastC._buildOwnData(t, Padding)

    #init hook necessaire pour info omp
    tmp   = Internal.getNodeFromName1(t, '.Solver#ownData')
    dtlocPy = Internal.getNodeFromName1(tmp, '.Solver#dtloc')  # noeud
    dtloc   = Internal.getValue(dtlocPy)                       # tab numpy
    nssiter = int(dtloc[0])

    zones = Internal.getZones(t)
    f_it = FastC.FIRST_IT
    if FastC.HOOK is None: FastC.HOOK = FastC.createWorkArrays__(zones, dtloc, f_it ); FastC.FIRST_IT = f_it

    # allocation d espace dans param_int pour stockage info openmp
    FastC._build_omp(t)

    zones_unstr =[]; zones_str = []; zones_ns = []; zones_lbm = []
    for z in zones:
        ztype = Internal.getValue(Internal.getNodeFromName(z, 'ZoneType'))
        if ztype == 'Structured':
            zones_str.append(z)
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
            if param_int[VSHARE.IFLOW] == 4:   #IFLOW=4
                zones_lbm.append(z)
            else:
                zones_ns.append(z)
        else: zones_unstr.append(z)

    # alloue metric: tijk, ventijk, ssiter_loc
    # init         : ssiter_loc
    metrics = allocate_metric(t, nghost)
    #print("apres alocate ", flush=True)

    metrics_str=[]; metrics_unstr=[]; metrics_ns=[]; metrics_lbm=[]
    c = 0
    for z in zones:
        ztype = Internal.getValue(Internal.getNodeFromName(z, 'ZoneType'))
        if ztype == 'Structured':
            metrics_str.append(metrics[c])
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]
            if param_int[VSHARE.IFLOW] == 4:   #IFLOW=4
                metrics_lbm.append(z)
            else:
                metrics_ns.append(z)
        else:
            metrics_unstr.append(metrics[c])
        c += 1
    # Contruction BC_int et BC_real pour CLa
    FastC._BCcompact(zones_str)

    #FastC._BCcompactNG(zones_unstr)

    #determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: range(22,1000)
    #print 'int(dtloc[0])= ', int(dtloc[0])
    if len(zones_str) != 0:
        for nstep in range(1, int(dtloc[0])+1):
            hook1       = FastC.HOOK.copy()
            hook1.update(FastC.fastc.souszones_list(zones_str, metrics_str, FastC.HOOK, 1, nstep, verbose) )

    #init metric
    FastC.fastc.init_metric(zones_str  , metrics_str  , hook1)
    #FastP.fastp.init_metric(zones_unstr, metrics_unstr, FastC.HOOK)

    # compact + align + init numa
    rmConsVars=True
    adjoint=Adjoint

    #print("avant primVavn resolved ", flush=True)
    t, FastC.FIRST_IT, zones2compact = FastC.createPrimVars(t, ompmode, rmConsVars, adjoint, lbmAJ=False)
    FastC.HOOK['FIRST_IT']= FastC.FIRST_IT

    if len(zones_str) == 0: hook1 = FastC.HOOK.copy(); hook1.update( {'lexit_lu':0, 'lssiter_verif':0})

    #compactage des champs en fonction option de calcul
    count = -1
    if ompmode == 1: count = 0
    for data in zones2compact:
        if ompmode == 1: count += 1
        zone    = data[0]
        varnames= data[1]
        for fields in varnames:
            #print("compact", fields)
            FastC._compact(zone, fields=fields, mode=count, dtloc=dtlocPy)
    for z in zones_unstr:
        _compact(zone, fields=None, mode=count, ParentElements=True)

    #print("apres compact", flush=True)
    #on recupere les zones a nouveau car create primvar rend zones caduc
    zones = Internal.getZones(t)
    #zones_str, zones_unstr, metrics_str, metrics_unstr = tri_zones( zones, metrics)
    infos_zones = tri_zones( zones, metrics)  #info-zone: dico 4 item contenant la list [ zones, metrics] pour Struct, unstruc, LBM, NS

    # init Q variables from macro variable if lbm et creation noeud pipeau density_P1 pour transfert
    flag_initprecise = 1
    if FASTLBM:
        FastLBM._init_Q(t, FastC.HOOK, flag_initprecise)

    #print("apres initQ", flush=True)
    #allocate tab ssor en structure
    zones_str   = infos_zones["struct"][0]
    metrics_str = infos_zones["struct"][1]
    ssors = FastS.fasts.allocate_ssor(zones_str, metrics_str, hook1)

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
    if ale and infos_ale is not None:
        print("ale actif. Teta et tetap=", infos_ale)
        teta = infos_ale[0];  tetap = infos_ale[1]
        FastC._motionlaw(t, teta, tetap)
        if len(zones_str) != 0: FastS.fasts.computePT_velocity_ale(zones_str,  metrics_str, FastC.HOOK, ompmode)
        #if len(zones_unstr) != 0: print("coder fastp.computePT_velocity_ale(zones_unstr,  metrics_unstr")
    #
    # Compactage arbre transfert et mise a jour FastC.HOOK
    #
    if tc is not None:
        X.miseAPlatDonorTree__(zones, tc, graph=graph, list_graph=list_graph)

        FastC.HOOK['param_int_tc'] = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
        param_real_tc = Internal.getNodeFromName1( tc, 'Parameter_real')
        if param_real_tc is not None: FastC.HOOK['param_real_tc']= param_real_tc[1]


        # Add linelets arrays in the FastC.HOOK for ODE-based WallModel (IBC)
        nbpts_linelets = 45
        #FastS._createTBLESA(tc,nbpts_linelets)
        FastS._createTBLESA(tc, h0=0.01, hn=-1, nbpts_linelets=nbpts_linelets)
    else:
        FastC.HOOK['param_real_tc'] = None
        FastC.HOOK['param_int_tc']  = None


    if ssors is not []:
        FastC.HOOK['ssors'] = ssors
    else:
        FastC.HOOK['ssors'] = None


    # Compactage arbre moyennes stat
    if tmy is not None:
        sol = Internal.getNodesFromName3(tmy, 'FlowSolution#Centers')
        var = Internal.getNodesFromType1(sol[0], 'DataArray_t')
        varmy = []
        for v in var: varmy.append('centers:'+v[0])
        FastC._compact(tmy, fields=varmy)

    #
    # remplissage ghostcells
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    timelevel_target = int(dtloc[4])

    #Determine variables a trasferer
    if len( infos_zones["LBM"][0]) ==0:
        # CALCUL FULL NS ==> on ne transfere que les variables macros
        vars=FastC.varsN
    elif len( infos_zones["NS"][0]) ==0 and len( infos_zones["unstruct"][0]) ==0:
        #CALCUL FULL LBM ==> on transfere les macros et les distributions
        vars=FastC.varsN + FastC.varsMLBM
    else:
        #Couplage NSLBM : on transfere les macro, les distributions et les gradients
        print("\n#------------Simulation couplee NS-LBM------------#")
        print('Nombre de zones LBM:',len(infos_zones['LBM'][0]))
        print('Nombre de zones NS:',len(infos_zones['NS'][0]))
        print('')
        vars = FastC.varsN + FastC.varsMLBM + FastC.varsS

    nstep            = 0
    nitrun           = 0

    # creation/init data dans param_int/real  pour interp temporelle LBM/NS et LBM/LBM
    init_interp = FastC._InterpTemporelcompact(t,tc)
    if FASTLBM:
        if not init_interp: FastLBM.fastaslbm.interp_temporel(zones, 1,1,-10)

    # data pour interp temporel couplage NS/LBM
    hook1 = FastC.HOOK.copy(); hook1.update( {'lexit_lu':0, 'lssiter_verif':0})

    #overset A basculer dans l'arbre pour eviter test systematique donnee statique
    overset = 0
    for z in zones:
        param_int = Internal.getNodeFromName(z, 'Parameter_int')[1]
        overset = max(overset,param_int[VSHARE.LBM_OVERSET])


    _fillGhostcells(zones, tc, infos_zones, timelevel_target, vars, nstep,  ompmode, hook1, overset=overset)
    if tc is not None: C._rmVars(tc, 'FlowSolution')

    #
    # initialisation Mut
    #
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    FastS.fasts._computePT_mut(infos_zones['struct'][0], infos_zones['struct'][1], hook1)

    return (t, tc, metrics)

#==============================================================================
# alloue retourne la metrique
#==============================================================================
def allocate_metric(t, nghost):
    zones        = Internal.getZones(t)
    dtloc        = Internal.getNodeFromName2(t, '.Solver#dtloc')
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

        ztype = Internal.getValue(Internal.getNodeFromName(z, 'ZoneType'))
        if ztype == 'Structured':
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')[1]

            if param_int[27] == 4:   #IFLOW=4
                if FASTLBM:
                    metrics.append(FastLBM.fastaslbm.allocate_metric(z, nssiter))
                else: raise ValueError("FastLBM needed but not installed.")
            else:
                metrics.append(FastS.fasts.allocate_metric(z, nssiter))
        #else: metrics.append(FastP.fastp.allocate_metric(z, nssiter))
    return metrics

#==============================================================================
def _compact(t, containers=[Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__], dtloc=None, fields=None, mode=None, init=True, ParentElements=False):
    if  mode is not None:
        if mode == -1: thread_numa = 0
        else: thread_numa = 1
    else: thread_numa = 0

    if dtloc is None:
        own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
        dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')    # noeud

    zones = Internal.getZones(t)
    ndom = 0
    for z in zones:
        if ParentElements:
            E_NG   = Internal.getNodeFromName(z   ,'NGonElements')
            FaceSch= Internal.getNodeFromName(z   ,'FaceScheduler')
            CellSch= Internal.getNodeFromName(z   ,'CellScheduler')
            a      = Internal.getNodeFromName(E_NG,'ParentElements')
            a1     = a[1]
            sh = a1.shape ; sh1= FaceSch[1].shape; sh2 = CellSch[1].shape

            eq = numpy.empty(sh , dtype=Internal.E_NpyInt, order='F')
            eq1= numpy.empty(sh1, dtype=Internal.E_NpyInt, order='F')
            eq2= numpy.empty(sh2, dtype=Internal.E_NpyInt, order='F')

            c = 0
            param_int = Internal.getNodeFromName2(z, 'Parameter_int')  # noeud
            if param_int is None:
                raise ValueError("_compact: Parameter_int is missing for zone %s."%z[0])
            opt = 1
            #if init: FastP.fastp.initNuma(a1, eq, eq1, eq2, param_int, c, thread_numa, opt)
            # Replace numpys with slice
            a[1]       = eq
            FaceSch[1] = eq1
            CellSch[1] = eq2

        else:
            opt = 0
            ztype = Internal.getValue(Internal.getNodeFromName(z, 'ZoneType'))
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
                        if ztype == 'Structured':
                            FastC.fastc.initNuma(ptr, eq,  param_int, dtloc, c, thread_numa, ndom)
                        #else: FastP.fastp.initNuma(ptr, eq, eq, eq, param_int, c, thread_numa, opt)

                    # Replace numpys with slice
                    a[1] = eq[c*(size)+c*param_int[1][66]:(c+1)*(size)+c*param_int[1][66]]
                    a[1] = a[1].reshape(sh, order='F')

                    c += 1
        ndom +=1

    return None

#==============================================================================
def _fillGhostcells(zones, tc, infos_zones, timelevel_target, vars, nstep, ompmode, hook1, nitmax=1, rk=1, exploc=0, num_passage=1, overset=0):

    #import KCore.test as test
    #test.testT(zones,  10)
    #test.testT(tc,  11)

    if hook1['lexit_lu'] == 0:

        #transfert
        if tc is not None :
            tc_compact = Internal.getNodeFromName1( tc, 'Parameter_real')
            #Si param_real n'existe pas, alors pas de raccord dans tc
            if tc_compact is not  None:

                couplageNSLBM = 0
                param_real    = tc_compact[1]
                param_int     = Internal.getNodeFromName1(tc, 'Parameter_int' )[1]
                zonesD        = Internal.getZones(tc)

                #-----------------------------------------------------------------
                #Definition de varType pour la fonction qui realise les transfert
                # = Cas Navier-Stokes
                if infos_zones["LBM"]==[[],[]]:
                    if hook1["neq_max"] == 5: varType =  2
                    else                    : varType = 21
                # = Cas LBM
                elif infos_zones["NS"]==[[],[]]:
                    if   hook1["neq_max"] == 19 and overset == 0: varType = 4 #LBM
                    elif hook1["neq_max"] == 19 and overset == 1:
                        varType = 41 #LBM avec Overset
                else: varType = 5; couplageNSLBM = 1 #Couplage NSLBM

                dtloc = hook1['dtloc']

                #-----------------------------------------------------------------
                #Copie dans le tc
                for v in vars:
                    if v in ['Q1','Q1_M1','corr_xx']: #Distributions => seules les zones LBM en ont
                        for l,z in enumerate(zones):
                            #print("cpvars",z[0])
                            if (z in infos_zones["LBM"][0]):
                                #print("addvars",z[0], v)
                                C._cpVars(z,'centers:'+v,zonesD[l],v)
                    else: C._cpVars(zones,'centers:'+v,zonesD,v)

                type_transfert = 2  # 0= ID uniquement, 1= IBC uniquement, 2= All
                no_transfert   = 1  # dans la list des transfert point a point
                isWireModel_int= 0
                #print("var", vars, varType)

                #t3 = C.newPyTree(['Base', zones])
                #t2 = Internal.copyTree(t3)
                #test.testT(t2, 200+nstep)

                nstep_loc = nstep
                if nstep==0: nstep_loc=1
                X.connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, timelevel_target, varType, type_transfert, no_transfert, nstep_loc, nitmax, rk, exploc, num_passage, isWireModel_int)

                #if couplageNSLBM==1 and nstep !=0 : fast.recuplbmns_(zones, zonesD, param_int, param_real, hook1, nitmax, 0, nstep_loc, ompmode, 1, 0)

                #t3 = C.newPyTree(['Base', zones])
                #t2 = Internal.copyTree(t3)
                #test.testT(t2, 700+nstep)

        #apply BC
        if exploc != 1:
            _applyBC(infos_zones, hook1, nstep, nitmax, var= vars)

        #t3 = C.newPyTree(['Base', zones])
        #t2 = Internal.copyTree(t3)
        #test.testT(t2, 800+nstep)

    return None

#==============================================================================
# compute in place
# graph is a dummy argument to be compatible with mpi version
#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None, layer="c", NIT=1):
    """Compute a given number of iterations."""

    import KCore.test as test

    own   = Internal.getNodeFromName1(t, '.Solver#ownData')  # noeud
    dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')    # noeud

    node = Internal.getNodeFromName1(t, '.Solver#define')
    omp_node = Internal.getNodeFromName1(node, 'omp_mode')
    ompmode  = OMP_MODE
    if  omp_node is not None: ompmode = Internal.getValue(omp_node)

    dtloc   = Internal.getValue(dtloc) # tab numpy
    nitmax  = int(dtloc[0])

    zones       = Internal.getZones(t)
    infos_zones = tri_zones( zones, metrics)

    #overset A basculer dans l'arbre pour eviter test systematique donnee statique
    overset = 0
    for z in zones:
        param_int = Internal.getNodeFromName(z, 'Parameter_int')[1]
        overset = max(overset,param_int[VSHARE.LBM_OVERSET])

    #if (model == 'CouplageNSLBM') and (len(infos_zones["LBM"][0])!=0) and (len(infos_zones["NS"][0])!=0):
    if len(infos_zones["LBM"][0])!=0 and len(infos_zones["NS"][0])!=0: flag_NSLBM = 1
    else: flag_NSLBM = 0

    #### a blinder...
    param_int_firstZone = Internal.getNodeFromName2( zones[0], 'Parameter_int' )[1]
    itypcp = param_int_firstZone[VSHARE.ITYPCP]
    exploc = param_int_firstZone[54]
    #### a blinder...

    if nitrun == 1: print('Info: using layer trans=%s (ompmode=%d)'%(layer, ompmode))

    if layer == "Python":

        if (exploc==1 or flag_NSLBM==1) and tc is not None:
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
            hook1.update(  FastC.fastc.souszones_list(zones, metrics, FastC.HOOK, nitrun, nstep, 0) )

            skip      = 0
            if hook1["lssiter_verif"] == 0 and nstep == nitmax and itypcp ==1: skip = 1

            # calcul Navier_stokes + appli CL
            if skip == 0:
                nstep_deb = nstep
                nstep_fin = nstep
                layer_mode= 0
                nit_c     = 1

                #if nstep==1 and flag_NSLBM==1:
                #   fast.interplbmns_(zones,zones_tc,param_int_tc,param_real_tc,hook1,nitmax,nitrun,0,ompmode,1,dest)

                #t3 = C.newPyTree(['Base', zones])
                #import KCore.test as test
                #t2 = Internal.copyTree(t)
                #test.testT(t2, 1100+nstep)

                param_int = Internal.getNodeFromName2(zones[0], 'Parameter_int' )[1]

                #t0=Time.time()
                sz= numpy.size(param_int)

                fast._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c, flag_NSLBM, hook1)
                #FastLBM.fastaslbm._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c, hook1)

                #if nstep ==3:
                #   C.convertPyTree2File(t,'2th_step3.cgns')
                #t2 = Internal.copyTree(t)
                #test.testT(t2, 1200+nstep)

                #t3 = C.newPyTree(['Base', zones])
                #t2 = Internal.copyTree(t3)
                #test.testT(t2, 1200+nstep)
                #if nstep==1: C.convertPyTree2File(t2,'verif.cgns')
                #import sys; sys.exit()


                #if flag_NSLBM==1: fast.interplbmns_(zones,zones_tc,param_int_tc,param_real_tc,hook1,nitmax,nitrun,nstep,ompmode,1,dest)
                #if flag_NSLBM==1: fast.interplbmns_(zones,zones_tc,param_int_tc,param_real_tc,hook1,nitmax,nitrun,nstep,1,dest)

                varsNS = FastC.varsP
                if nstep%2 == 0 and itypcp == 2 : varsNS = FastC.varsN  # Choix du tableau pour application transfer et BC

                #Definition de vars pour la fonction qui realise les transfert
                if   len(infos_zones['LBM'][0])==0: vars = varsNS                         # full NS
                elif len(infos_zones['NS'][0] )==0:                                        # full LBM
                    vars = FastC.varsN + FastC.varsPLBM                                   # LBM match
                    if hook1["neq_max"] == 19 and overset==1 :                            # LBM avec Overset
                        vars = FastC.varsN + FastC.varsPLBM + FastC.varsS + FastC.varsPSI
                else:
                    vars = varsNS + FastC.varsPLBM + FastC.varsS

                timelevel_target = int(dtloc[4])

                #t3 = C.newPyTree(['Base', zones])
                #t2 = Internal.copyTree(t3)
                #test.testT(t2, 1300+nstep)

                _fillGhostcells(zones, tc,  infos_zones, timelevel_target, vars, nstep, ompmode, hook1, overset=overset, nitmax=nitmax)

                #t3 = C.newPyTree(['Base', zones])
                #t2 = Internal.copyTree(t3)
                #test.testT(t2, 1400+nstep)

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

        fast._computePT(zones, metrics, nitrun, nstep_deb, nstep_fin, layer_mode, nit_c , flag_NSLBM,  FastC.HOOK)



    zones_ns    = infos_zones["NS"][0]
    metrics_ns  = infos_zones["NS"][1]
    zones_lbm   = infos_zones["LBM"][0]
    metrics_lbm = infos_zones["LBM"][1]

    #switch pointer a la fin du pas de temps
    if exploc == 1 and tc is not None:
        if layer == 'Python':
            FastC.switchPointers__(zones_ns, 1, 3)
        else:
            FastC.switchPointers3__(zones_ns,nitmax)
    else:
        case = NIT%3
        if case != 0 and itypcp < 2: FastC.switchPointers__(zones_ns, case)
        if case != 0 and itypcp ==2: FastC.switchPointers__(zones_ns, case, nitmax%2+2)

    if FASTLBM:
        if  NIT%2 != 0:
            FastC.switchPointersLBM__(zones_lbm, FastLBM.NQ, dtloc)

    maxlevel   = dtloc[9]
    nitCyclLBM = 2**(maxlevel-1)
    dtloc[10] += 1    #itCycl_lbm
    if dtloc[10] == nitCyclLBM: dtloc[10]=0


    # flag pour derivee temporelle 1er pas de temps implicit
    FastC.HOOK["FIRST_IT"]  = 1
    FastC.FIRST_IT          = 1

    return None

#==============================================================================
# applyBC
#==============================================================================
def _applyBC(infos_zones, hook1, nstep, nitmax, var=["Density","Q1"]):

    if nstep==0: nstep_NS = 1; nstep_LBM = 1
    else       : nstep_NS = nstep; nstep_LBM = nstep

    nstep_bc = 0
    if nitmax==3: nstep_bc = 3
    elif nitmax > 3: nstep_bc = nitmax-1

    #Variables on which the BCs are applied
    varns = var[0]; varlbm = "Q1_M1"
    if len(var)>=2: varlbm = var[1]
    if len(var)==4: varlbm = var[-1]

    #print("fast: varBC NS-LBM", varns, varlbm)
    #On applique les BC selon les domaines en presence
    #Navier Stokes structure : FastS
    if len(infos_zones['NS_str'][0]) != 0:
        FastS.fasts._applyBC(infos_zones["NS_str"][0]  , infos_zones["NS_str"][1]  , hook1, nstep_NS,  varns  )
        #print('FastS BC applied')
    #LBM : FastLBM
    if FASTLBM:
        if len(infos_zones['LBM'][0]) != 0 and (nstep==0 or nstep==nstep_bc):
            FastLBM.fastaslbm._applyBC(infos_zones["LBM"][0] , infos_zones["LBM"][1]  , hook1, nstep_LBM, varlbm )
            #if nstep==0: print('FastLBM BC applied 1')
            #else: print('FastLBM BC applied 2')

    #Navier Stokes non structure : FastP
    #if len(infos_zones['NS_ustr'][0]) != 0:
    #   FastP.fastp._applyBC(infos_zones["unstruct"][0], infos_zones["unstruct"][1]  , hook1, nstep_NS, ompmode, varns  )

    return None

#==============================================================================
# renvoi list zone structure et polyhedrique a partir arbre hybride
#==============================================================================
def tri_zones(zones, metrics):

    infos_zones = {}
    zones_unstr =[]; zones_str =[]; zones_ns =[]; zones_lbm =[]; zones_ns_str =[]
    metrics_unstr =[]; metrics_str =[]; metrics_ns =[]; metrics_lbm =[]; metrics_ns_str = []
    ind_LBM = []

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
                ind_LBM.append(c)
            else:
                zones_ns.append(z)
                zones_ns_str.append(z)
                metrics_ns.append( metrics[c] )
                metrics_ns_str.append( metrics[c])
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
    infos_zones['NS_str'  ]=[zones_ns_str, metrics_ns_str]
    infos_zones['NS_ustr' ]=[zones_unstr, metrics_unstr]
    infos_zones['iLBM'    ]=ind_LBM
    return infos_zones

#==============================================================================
# Display
#==============================================================================
def display_temporal_criteria(t, metrics, nitrun, format=None, gmres=None,verbose='firstlast'):
    zones        = Internal.getZones(t)
    dtloc        = Internal.getNodeFromName2(t, '.Solver#dtloc')
    dtloc_numpy  = Internal.getValue(dtloc)
    nssiter      = int(dtloc_numpy[0])
    nzones	 = len(zones)

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

    infos_zones = tri_zones( zones, metrics)  #info-zone: dico 4 keys contenant la list [ zones, metrics] pour Struct, unstruc, LBM, NS

    iverb = 0
    if verbose != 'firstlast': iverb=2

    residu = FastS.fasts.display_ss_iteration(infos_zones['struct'][0], infos_zones['struct'][1], cvg_numpy, nitrun, nssiter, lft, iverb)
    #residu = FastP.fastp.display_ss_iteration( infos_zones['unstruct'][0], infos_zones['unstruct'][1], cvg_numpy, nitrun, nssiter, lft)

    if gmres is None: return None
    else: return residu


#==============================================================================
# Initialisation parametres interpolation LBM --> NS : passage en revue des
# raccord LBM -> NS, mise en place stockage et donnees pour l'interp
#==============================================================================
def prepare_interpolation(t,tc):
    zones         = Internal.getZones(t)
    zones_tc      = Internal.getZones(tc)
    nzones        = len(zones)

    neq = 5 #ou 6 pour RANS ?
    timelevels = 3
    nb_gc = 2

    # 0 : nb_rac, 1 : taille_rac_max, 2-37 : range_save, 38-73 : gc_fill
    dim_interp = 6*6 + 6*6 + 1 + 1

    rac_insta = {}
    size = 0
    print('===================Prepare interpolation=================')
    for z in zones:
        #print("Zone =",z[0])
        dim = Internal.getZoneDim(z); #print(dim)
        model="NSLaminar"
        a = Internal.getNodeFromName2(z,'GoverningEquations')
        if a is not None : model = Internal.getValue(a)
        nb_interp = 0

        if model  == 'LBMLaminar':
            rac_insta[z[0]] = []
            rac_insta[z[0]].append(nb_interp)
            rac_insta[z[0]].append(nb_interp)
            for bc in Internal.getNodesFromType2(z,'BC_t'):
                v = Internal.getValue(bc); #print(v)
                if v=='BCReconsLBM':
                    nb_interp += 1
                    node = Internal.getNodeFromName(bc,'PointRange')
                    node = sum(Internal.getValue(node).tolist(),[])
                    pt_to_save, gc_to_fill = get_InterpDataFromNode(node,dim,nb_gc,'BCReconsLBM')
                    rac_insta[z[0]].append(pt_to_save); rac_insta[z[0]].append(gc_to_fill)
                    taille_rac = (pt_to_save[1]-pt_to_save[0]+1)*(pt_to_save[3]-pt_to_save[2]+1)*(pt_to_save[5]-pt_to_save[4]+1)
                    #print(taille_rac)
                    size = max(size,taille_rac)
                # elif v=='BCadimcoins':
                #    nb_interp += 1
                #    node = Internal.getNodeFromName(bc,'PointRange')
                #    node = sum(Internal.getValue(node).tolist(),[])
                #    pt_to_save, gc_to_fill = get_InterpDataFromNode(node,dim,nb_gc,'BCadimcoins')
                #    rac_insta[z[0]].append(pt_to_save); rac_insta[z[0]].append(gc_to_fill)
                #    taille_rac = (pt_to_save[1]-pt_to_save[0]+1)*(pt_to_save[3]-pt_to_save[2]+1)*(pt_to_save[5]-pt_to_save[4]+1)
                #    size = max(size,taille_rac)
            rac_insta[z[0]][0] = nb_interp; #print(nb_interp)
    for z in zones:
        if z[0] in rac_insta.keys():
            rac_insta[z[0]][1] = size
    print(rac_insta)

    nb_rac_tot = 0
    size_max = size
    for key in rac_insta:
        info_rac = rac_insta[key]
        nb_rac_tot += info_rac[0]
        interp_info = numpy.zeros(dim_interp, dtype=Internal.E_NpyInt)
        interp_info[0] = info_rac[0] #Nombre de raccors LBM -> NS pour la zone
        interp_info[1] = size_max*neq*timelevels
        for i in range(info_rac[0]):
            shift_pt_rac = 2 + i*2
            shift_gc_rac = 2 + i*2 + 1
            shift_pt_tab = 2 + i*6
            shift_gc_tab = 2 + 6*6 + i*6
            for j in range(6):
                interp_info[shift_pt_tab+j] = info_rac[shift_pt_rac][j]
                interp_info[shift_gc_tab+j] = info_rac[shift_gc_rac][j]
        #print(interp_info)
        z = Internal.getNodeFromName(t,key)
        o = Internal.createUniqueChild(z,'.Solver#ownData','UserDefinedData_t')
        Internal.createUniqueChild(o,'Interp_data','DataArray_t',interp_info)

    # Tableau de stockage pour interpolation LBM->NS
    dim_tabinterp = 5*3*nb_rac_tot*size_max
    print('Dim Tab interp NS/LBM temporel:',dim_tabinterp)
    tab_interp = numpy.empty(dim_tabinterp, dtype=numpy.float64)
    FastC.HOOK['tab_interp'] = tab_interp

    return None

#==============================================================================
def get_InterpDataFromNode(node,dim,nb_gc,bc_name):
    #On calcule les range des frontieres de la zone
    imin = [     1,      1,      1, dim[2],      1, dim[3]]
    imax = [dim[1], dim[1],      1, dim[2],      1, dim[3]]
    jmin = [     1, dim[1],      1,      1,      1, dim[3]]
    jmax = [     1, dim[1], dim[2], dim[2],      1, dim[3]]
    kmin = [     1, dim[1],      1, dim[2],      1,      1]
    kmax = [     1, dim[1],      1, dim[2], dim[3], dim[3]]
    #print(node)
    #print(imin, imax, jmin, jmax, kmin, kmax)
    if node == imin and bc_name=='BCReconsLBM':
        pt_to_save = [1,1+(nb_gc-1),1,(dim[2]-1)-2*nb_gc,1,(dim[3]-1)-2*nb_gc]
        gc_to_fill = [pt_to_save[0]-nb_gc,pt_to_save[1]-nb_gc,pt_to_save[2]-nb_gc,pt_to_save[3]+nb_gc,pt_to_save[4]-nb_gc,pt_to_save[5]+nb_gc]
        #print(pt_to_save, gc_to_fill)
    elif node == imin and bc_name=='BCadimcoins':
        #print("Node = imin")
        pt_to_save = [0,0,0,0,0,0]
        gc_to_fill = [-1,0,-1,0,1-nb_gc,(dim[3]-1)-nb_gc]
        #print(gc_to_fill)
    elif node == imax and bc_name=='BCReconsLBM':
        pt_to_save = [(dim[1]-1)-2*nb_gc-1,(dim[1]-1)-2*nb_gc,1,(dim[2]-1)-2*nb_gc,1,(dim[3]-1)-2*nb_gc]
        gc_to_fill = [pt_to_save[0]+nb_gc,pt_to_save[1]+nb_gc,pt_to_save[2]-nb_gc,pt_to_save[3]+nb_gc,pt_to_save[4]-nb_gc,pt_to_save[5]+nb_gc]
        #print(pt_to_save,gc_to_fill)
    elif node == imax and bc_name=='BCadimcoins':
        #print("Node = imax")
        pt_to_save = [0,0,0,0,0,0]
        gc_to_fill = [(dim[1]-1)-nb_gc-1,(dim[1]-1)-nb_gc,(dim[2]-1)-nb_gc-1,(dim[2]-1)-nb_gc,1-nb_gc,(dim[3]-1)-nb_gc]
    elif node == jmin and bc_name=='BCReconsLBM':
        pt_to_save = [1,(dim[1]-1)-2*nb_gc,1,1+(nb_gc-1),1,(dim[3]-1)-2*nb_gc]
        gc_to_fill = [pt_to_save[0]-nb_gc,pt_to_save[1]+nb_gc,pt_to_save[2]-nb_gc,pt_to_save[3]-nb_gc,pt_to_save[4]-nb_gc,pt_to_save[5]+nb_gc]
        #print(pt_to_save,gc_to_fill)
    elif node == jmin and bc_name=='BCadimcoins':
        #print("Node = jmin")
        pt_to_save = [0,0,0,0,0,0]
        gc_to_fill = [(dim[1]-1)-nb_gc-1,(dim[1]-1)-nb_gc,-1,0,1-nb_gc,(dim[3]-1)-nb_gc]
        #print(gc_to_fill)
    elif node == jmax and bc_name=='BCReconsLBM':
        pt_to_save = [1,(dim[1]-1)-2*nb_gc,(dim[2]-1)-2*nb_gc-1,(dim[2]-1)-2*nb_gc,1,(dim[3]-1)-2*nb_gc]
        gc_to_fill = [pt_to_save[0]-nb_gc,pt_to_save[1]+nb_gc,pt_to_save[2]+nb_gc,pt_to_save[3]+nb_gc,pt_to_save[4]-nb_gc,pt_to_save[5]+nb_gc]
        #print(pt_to_save,gc_to_fill)
    elif node == jmax and bc_name=='BCadimcoins':
        #print("Node = jmax")
        pt_to_save = [0,0,0,0,0,0]
        gc_to_fill = [-1,0,(dim[2]-1)-nb_gc-1,(dim[2]-1)-nb_gc,1-nb_gc,(dim[3]-1)-nb_gc]
        #print(gc_to_fill)
    elif node == kmin and bc_name=='BCReconsLBM':
        pt_to_save = [1,(dim[1]-1)-2*nb_gc,1,(dim[2]-1)-2*nb_gc,1,1+(nb_gc-1)]
        gc_to_fill = [pt_to_save[0]-nb_gc,pt_to_save[1]+nb_gc,pt_to_save[2]-nb_gc,pt_to_save[3]+nb_gc,pt_to_save[4]-nb_gc,pt_to_save[5]-nb_gc]
    elif node == kmax and bc_name=='BCReconsLBM':
        pt_to_save = [1,(dim[1]-1)-2*nb_gc,1,(dim[2]-1)-2*nb_gc,(dim[3]-1)-2*nb_gc-1,(dim[3]-1)-2*nb_gc]
        gc_to_fill = [pt_to_save[0]-nb_gc,pt_to_save[1]+nb_gc,pt_to_save[2]-nb_gc,pt_to_save[3]-nb_gc,pt_to_save[4]+nb_gc,pt_to_save[5]+nb_gc]
    else:
        print('Unknown boundary for interpolation')

    return pt_to_save, gc_to_fill

#==============================================================================
def recup_fgc(zones):
    for z in zones:
        paramint = Internal.getNodeFromName2(z, 'Parameter_int')[1]
        if paramint[27]==4:
            C._initVars(z,'{centers:Density}={centers:Density_P1}')
            C._initVars(z,'{centers:VelocityX}={centers:VelocityX_P1}')
            C._initVars(z,'{centers:VelocityY}={centers:VelocityY_P1}')
            C._initVars(z,'{centers:VelocityZ}={centers:VelocityZ_P1}')
