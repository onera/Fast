# FastS + MPI
import PyTree
import fasts
from PyTree import display_temporal_criteria, createConvergenceHistory, extractConvergenceHistory, createStressNodes, _computeStress, metric, createPrimVars, _createPrimVars, createStatNodes, _computeStats, initStats, _computeEnstrophy, _computeVariables, _computeGrad, _compact, _applyBC, _buildOwnData,  metric, _BCcompact, _computeVelocityAle, checkBalance, itt

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Connector.Mpi as Xmpi
    import Connector.PyTree as X
    import Connector
    import Fast.Internal as FastI
except:
    raise ImportError("FastS: requires Converter and Connector modules.")

#==============================================================================
def _compute(t, metrics, nitrun, tc=None, graph=None):
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

    for nstep in xrange(1, nitmax+1): # pas RK ou ssiterations
        # determination taille des zones a integrer (implicit ou explicit local)
        hook1 = PyTree.HOOK + fasts.souszones_list(zones, metrics, PyTree.HOOK, nitrun, nstep)
        nidom_loc = hook1[11]

        # hook1[10] = nombre equations
        # hook1[11] = nidom_lu
        # hook1[12] = lskip_lu
        # hook1[13] = lssiter_verif

        skip = 0
        if (hook1[13] == 0 and nstep == nitmax and itypcp ==1): skip = 1

        # calcul Navier Stokes + appli CL
        if nidom_loc > 0 and skip ==0:

            # Navier-Stokes
            fasts._computePT(zones, metrics, nitrun, nstep, omp_mode, hook1)

            # Ghostcell
            vars = PyTree.varsP
            if  nstep == 2 and itypcp == 2 : vars = PyTree.varsN  # Choix du tableau pour application transfer et BC
            _fillGhostcells(zones, tc, metrics, nitrun, vars, nstep, hook1, graphID, graphIBCD, procDict)

    # switch pointers
    FastI.switchPointers__(zones, orderRk)
    # flag pour derivee temporelle 1er pas de temps implicit
    PyTree.HOOK[9]  = 1
    PyTree.FIRST_IT = 1
    return None

#==============================================================================
def _fillGhostcells(zones, tc, metrics, nitrun, vars, nstep, hook1, graphID, graphIBCD, procDict): 
    
   # hook1[10] = nombre equations max
   # hook1[11] = nidom_lu
   # hook1[12] = lskip_lu
   # hook1[13] = lssiter_verif
   if hook1[12] ==0:

       #transfert
       if tc is not None :
           tc_compact = Internal.getNodeFromName1(tc, 'Parameter_real')
           if tc_compact is not None:

              param_real= tc_compact[1]
              param_int = Internal.getNodeFromName1( tc, 'Parameter_int' )[1]
              zonesD    = Internal.getZones(tc)

              if hook1[10] == 5: varType = 2
              else             : varType = 21

              bcType = PyTree.HOOKIBC[0]; Gamma=PyTree.HOOKIBC[1]; Cv=PyTree.HOOKIBC[2]; Mus=PyTree.HOOKIBC[3]; Cs=PyTree.HOOKIBC[4]; Ts=PyTree.HOOKIBC[5]

              if nstep <= 3: 
                 for v in vars: C._cpVars(zones, 'centers:'+v, zonesD, v)

              #recuperation Nb pas instationnnaire dans tc
              type_transfert = 1  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, nitrun,
                                        bcType=bcType, varType=varType, compact=1,
                                        graph=graphIBCD, procDict=procDict, 
                                        Gamma=Gamma, Cv=Cv, MuS=Mus, Cs=Cs, Ts=Ts)
              type_transfert = 0  # 0= ID uniquememnt, 1= IBC uniquememnt, 2= All 
              Xmpi.__setInterpTransfers(zones , zonesD, vars, param_int, param_real, type_transfert, nitrun,
                                        varType=varType, compact=1, graph=graphID, procDict=procDict)

       #apply BC
       _applyBC(zones, metrics, var=vars[0])

   return None
#==============================================================================
def warmup(t, tc, graph=None, infos_ale=None):
    #global FIRST_IT, HOOK, HOOKIBC

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
    # Calul de la metric: tijk, ventijk, ssiter_loc
    metrics = metric(t)
    # Contruction BC_int et BC_real pour CL
    _BCcompact(t) 
    # compact + align + init numa
    t = createPrimVars(t, omp_mode)

    # determination taille des zones a integrer (implicit ou explicit local)
    #evite probleme si boucle en temps ne commence pas a it=0 ou it=1. ex: xrange(22,1000)
    dtloc = Internal.getNodeFromName3(t, '.Solver#dtloc')  # noeud
    dtloc = Internal.getValue(dtloc)                       # tab numpy
    zones = Internal.getZones(t)

    #Allocation HOOK
    f_it = PyTree.FIRST_IT
    if PyTree.HOOK is None: PyTree.HOOK = FastI.createWorkArrays__(zones, dtloc, f_it ); PyTree.FIRST_IT = f_it
    for nstep in xrange(1, int(dtloc[0])+1): hook1 = PyTree.HOOK + fasts.souszones_list(zones, metrics, PyTree.HOOK, 1, nstep)

    #Allocation HOOKIBC
    if PyTree.HOOKIBC is None: PyTree.HOOKIBC = FastI.getIBCInfo__(t)

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
    if(ale == True and infos_ale is not None):
        print "ale actif. Teta et tetap=", infos_ale
        teta = infos_ale[0];  tetap = infos_ale[1]
        _motionlaw(t, teta, tetap)
        _computeVelocityAle(t,metrics)

    #
    # Compactage arbre transfert
    #
    if tc is not None:
       if graph is not None: 
          g = graph['procDict']
          l = graph['procList']
       else: 
          g = None; l = None

       X.miseAPlatDonnorTree__(zones, tc, procDict=g, procList=l)

    #
    # remplissage ghostcells
    #
    hook1[12] = 0
    nstep     = 1
    nitrun    = 0
    if infos_ale is not None and len(infos_ale) == 3: nitrun = infos_ale[2]
    _fillGhostcells(zones, tc, metrics, nitrun, ['Density'], nstep, hook1,graphID, graphIBCD, procDict) 
    
    return (t, tc, metrics)
