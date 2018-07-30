"""Common functions for FAST solvers.
"""
import Fast
import numpy
__version__ = Fast.__version__

try: import Fast.Internal as FastI
except: import Internal as FastI
from Internal import _setNum2Zones, _setNum2Base, setNum2Zones, setNum2Base
try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
except:
    raise ImportError("Fast.PyTree: requires Converter module.")

# Solveurs disponibles
FASTP = True
try: import FastP.PyTree as FastP
except: FASTP = False
FASTS = True
try: import FastS.PyTree as FastS
except: FASTS = False

#==============================================================================
# metric
# Retourne t, metrics
#==============================================================================
def metric(t):
    zones = Internal.getZones(t)
    metrics = []
    for z in zones:
        solver = getValueFromTag(z, 'solver')
        if solver == 'FastP' and FASTP: 
            #metric = FastP.metric(z); 
            metrics.append(None)
        elif solver == 'FastS' and FASTS: 
            #metric = FastS.metric(z); 
            #metrics.append(metric)
            metrics.append(None)
        else: metrics.append(None)
    return t, metrics

#==============================================================================
# compute
#==============================================================================
def compute(t, metrics):
    tp = Internal.copyRef(t)
    _compute(tp, metrics)
    return tp

def _compute(t, metrics):
    zones = Internal.getZones(t)
    i = 0
    for z in zones:
        solver = getValueFromTag(z, 'solver')
        metric = metrics[i]
        if solver == 'FastP' and FASTP: 
            #FastP._computeZone(z)
            pass
        elif solver == 'FastS' and FASTS: 
            #FastS._computeZone(z, metric)
            pass
        i += 1
    return None

#==============================================================================
# applyBC
#==============================================================================
def applyBC(t, metrics, topTree=None):
    tp = Internal.copyRef(t)
    _applyBC(tp, metrics, topTree)
    return tp

def _applyBC(t, metrics, topTree=None):
    if topTree is None: top = t
    else: top = topTree
    state = FastI.getRefState__(top)
    zones = Internal.getZones(t)
    i = 0
    for z in zones:
        solver = getValueFromTag(z, 'solver')
        metric = metrics[i]
        if solver == 'FastP' and FASTP: 
          #FastP._applyBCZone(z, state)
          pass
        elif solver == 'FastS' and FASTS: 
          #FastS._applyBCZone(z, metric, state)
          pass
        i += 1
    return None

#==============================================================================
# Calcul et retourne la metrique
#==============================================================================
#def metric(t):
#    #global FIRST_IT
#    zones        = Internal.getNodesFromType2(t, 'Zone_t')
#    dtloc        = Internal.getNodeFromName3(t, '.Solver#dtloc')
#    dtloc_numpy  = Internal.getValue(dtloc)
#    nssiter      = int(dtloc_numpy[0])
#    #FIRST_IT     = int(dtloc_numpy[2])
#    metrics=[]; motion ='none'
#    for z in zones:
#        b = Internal.getNodeFromName2(z, 'motion')
#        if b is not None: motion = Internal.getValue(b)
#        num = Internal.getNodeFromName1(z, '.Solver#ownData')
#        if num is None:
#            raise ValueError("metric: numerics is missing for zone %s."%z[0])
#        if motion == 'rigid':
#            grids = Internal.getNodesFromType1(z, 'GridCoordinates_t')
#            if len(grids) == 1:
#               grid_init = Internal.copyTree(grids[0])
#               grid_init[0] = 'GridCoordinates#Init'
#               Internal.addChild(z, grid_init, pos=1) # first
#        b    = Internal.getNodeFromName1(z, 'ZoneType')
#        topo = Internal.getValue(b)
#        #if(topo == 'Structured'): metrics.append(FastS.fasts.metric(z, nssiter))
#        #else                    : metrics.append(FastP.fastp.metric(z, nssiter))
#    return metrics


#==============================================================================
# Retourne un dictionnaire num a partir des donnees solveur d'une zone
#==============================================================================
def getNumFromTag(t):
    num = {}
    z = Internal.getNodeFromType2(t, 'Zone_t')
    if z is None: return num
    node = Internal.getNodeFromName1(z, '.Solver#define')
    if node is None: return num
    for n in node[2]: num[n[0]] = Internal.getValue(n)
    return num

#==============================================================================
# Retourne la valeur d'un champ de tag
# IN: name: champ voulu (solver...)
#==============================================================================
def getValueFromTag(t, name):
    z = Internal.getNodeFromType2(t, 'Zone_t')
    if z is None: return None
    node = Internal.getNodeFromName1(z, '.Solver#define')
    if node is None: return None
    node = Internal.getNodeFromName1(node, name)
    if node is None: return None
    val = Internal.getValue(node)
    return val

#==============================================================================
# Cette methode ne sert que pour le cas 2D (avec plusieurs plans)
# Elle fait passer le cellN des elements exterieurs a 4
# Elle n'agit que pour le solveur P
#==============================================================================
def adaptCellNForGC(t, depth, metrics):
    """Modify cellN for elements in the neighbourhood of depth layers of elements to the external borders."""
    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp)
    i = 0
    for z in zones:
        solver = getValueFromTag(z, 'solver')
        if solver == 'FastP' and FASTP == True: 
            FastP._adaptCellNForGCZone(z, depth, metrics[i])
        i += 1
    return tp

#==============================================================================
# Load t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
#==============================================================================
def load(fileName='t', fileNameC='tc', fileNameS='tstat', split='single',
         restart=False, NP=0):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension
    baseNameC = os.path.basename(fileNameC)
    baseNameC = os.path.splitext(baseNameC)[0] # name without extension
    fileNameC = os.path.splitext(fileNameC)[0] # full path without extension
    baseNameS = os.path.basename(fileNameS)
    baseNameS = os.path.splitext(baseNameS)[0] # name without extension
    fileNameS = os.path.splitext(fileNameS)[0] # full path without extension

    graph = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if NP > 0: # mpi run
        import Converter.Mpi as Cmpi
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            # Load connect (tc)
            FILE = fileNameC+'.cgns'
            tc = Cmpi.convertFile2SkeletonTree(FILE)
            tc = Cmpi.readZones(tc, FILE, rank=rank)
            graphID = Cmpi.computeGraph(tc, type='ID')
            graphIBCD = Cmpi.computeGraph(tc, type='IBCD')
            procDict = D2.getProcDict(tc)
            procList = D2.getProcList(tc, sort=True)
            graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
            tc = Cmpi.convert2PartialTree(tc, rank=rank)
            # Load data (t)
            FILE = fileName+'.cgns'
            if restart and os.access('restart.cgns', os.F_OK):
                FILE = 'restart.cgns'
            if os.access(FILE, os.F_OK):
                t = Cmpi.convertFile2SkeletonTree(FILE)
                t = Cmpi.readZones(t, FILE, rank=rank)
                t = Cmpi.convert2PartialTree(t)
            else: t = None
            # Load stat (ts)
            FILE = fileNameS+'.cgns'
            if os.access(FILE, os.F_OK):
                ts = Cmpi.convertFile2SkeletonTree(FILE)
                ts = Cmpi.readZones(ts, FILE, rank=rank)
                ts = Cmpi.convert2PartialTree(ts)
            else: ts = None

        else: # un fichier loade par proc
            # Try to load graph
            if os.access('%s/graph.pk'%fileName, os.R_OK):
                import cPickle as pickle
                file = open('%s/graph.pk'%fileName, 'rb')
                graph = pickle.load(file)
                file.close()
            FILE = '%s/%s_%d.cgns'%(fileNameC, baseNameC, rank)
            if os.access(FILE, os.F_OK): tc = C.convertFile2PyTree(FILE)
            else:
                FILE = fileNameC+'.cgns'
                tc = Cmpi.convertFile2SkeletonTree(FILE)
                tc = Cmpi.readZones(tc, FILE, rank=rank)        
                graphID = Cmpi.computeGraph(tc, type='ID',reduction=False)
                graphIBCD = Cmpi.computeGraph(tc, type='IBCD',reduction=False)
                procDict = D2.getProcDict(tc)
                procList = D2.getProcList(tc, sort=True)
                graph = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                tc = Cmpi.convert2PartialTree(tc)
            # Load data (t)
            FILE = '%s/%s_%d.cgns'%(fileName, baseName, rank)
            if restart and os.access('restart/restart_%d.cgns'%rank, os.F_OK):
                FILE = 'restart/restart_%d.cgns'%rank
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
            # Load stat (ts)
            FILE = '%s/%s_%d.cgns'%(fileNameS, baseNameS, rank)
            if os.access(FILE, os.F_OK): ts = C.convertFile2PyTree(FILE)
            else: ts = None
    else: # sequential run
        if split == 'single':
            # Load Connectivity (tc)
            FILE = fileNameC+'.cgns'
            if os.access(FILE, os.F_OK): tc = C.convertFile2PyTree(FILE)
            else: tc = None
            # Load Data (t)
            FILE = fileName+'.cgns'
            if restart and os.access('restart.cgns', os.F_OK):
                FILE = 'restart.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
            # Load Stat (ts)
            FILE = fileNameS+'.cgns'
            if os.access(FILE, os.F_OK): ts = C.convertFile2PyTree(FILE)
            else: ts = None
        else: # multiple
            # Load connectivity (tc)
            ret = 1; no = 0; tc = []
            while ret == 1:
                FILE = '%s/%s_%d.cgns'%(fileNameC, baseNameC, no)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: tc.append(C.convertFile2PyTree(FILE))
                no += 1
            if tc != []: tc = Internal.merge(tc)
            else: 
                FILE = fileNameC+'.cgns'
                if os.access(FILE, os.F_OK): tc = C.convertFile2PyTree(FILE)
                else: tc = None
            # Load data (t)
            ret = 1; no = 0; t = []
            while ret == 1:
                FILEIN = '%s/%s_%d.cgns'%(fileName, baseName, no)
                if os.access(FILEIN, os.F_OK): FILE = FILEIN
                elif restart and os.access('restart/restart_%d.cgns'%no, os.F_OK): 
                    FILE = 'restart/restart_%d.cgns'%no
                else: ret = 0
                if ret == 1: t.append(C.convertFile2PyTree(FILE))
                no += 1
            if t != []: t = Internal.merge(t)
            else: t = None
            # Load stat (ts)
            ret = 1; no = 0; ts = []
            while ret == 1:
                FILE = '%s/%s_%d.cgns'%(fileNameS, baseNameS, no)
                if not os.access(FILEIN, os.F_OK): ret = 0 
                if ret == 1: ts.append(C.convertFile2PyTree(FILE))
                no += 1
            if ts != []: ts = Internal.merge(ts)
            else: ts = None
    return t, tc, ts, graph

#==============================================================================
# save t
# IN: NP: 0 (seq run), >0 (mpi run, distributed), <0 (seq, save multiple)
# IN: split: 'single', 'multiple'
# IN: fileName: name of output file or dir
# single: write restart.cgns (full)
# multiple: write restart/restart_1.cgns, ... (partial trees)
#==============================================================================
def save(t, fileName='restart', split='single',
         temporal_scheme='implicit', NP=0):
    """Save tree and connectivity tree."""
    # Rip file ext if any
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    # Rip some useless data (FastS)
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')
    zones = Internal.getZones(t2)

    flowsol = Internal.getNodeFromName1(zones[0], 'FlowSolution#Centers')
    if flowsol is not None:
        vars    = Internal.getNodesFromType1(flowsol, 'DataArray_t')
        for var in vars:
            if (('Kry' in var[0]) and not('Kry_0' in var[0])):
                C._rmVars(t2, 'centers:'+var[0])

    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#define')
        if node is not None:
            node = Internal.getNodeFromName1(node, 'temporal_scheme')
            if node is not None:
                integ = Internal.getValue(node)
                if integ == 'explicit':
                    C._rmVars(z, 'centers:Density_M1')
                    C._rmVars(z, 'centers:VelocityX_M1')
                    C._rmVars(z, 'centers:VelocityY_M1')
                    C._rmVars(z, 'centers:VelocityZ_M1' )
                    C._rmVars(z, 'centers:Temperature_M1')
                    C._rmVars(z, 'centers:TurbulentSANuTilde_M1')

    # save
    if NP > 0: # mpi run
        import Converter.Mpi as Cmpi
        if split == 'single':  # output in a single file
            Cmpi.convertPyTree2File(t2, fileName+'.cgns')
        else:
            rank = Cmpi.rank
            C.convertPyTree2File(t2, '%s/%s_%d.cgns'%(fileName,baseName,rank))
            # Rebuild graph
            # skeleton -> gather -> merge -> graph
    else: # sequential run
        if split == 'single': 
            C.convertPyTree2File(t2, fileName+'.cgns')
        else:
            # Get and save graph
            import Converter.Mpi as Cmpi
            import Distributor2.PyTree as D2
            graphID = Cmpi.computeGraph(t2, type='ID')
            graphIBCD = Cmpi.computeGraph(t2, type='IBCD')
            procDict = D2.getProcDict(t2)
            procList = D2.getProcList(t2,  sort=True)
            objet = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList}

            # Rebuild local trees
            if not os.path.exists(fileName): os.makedirs(fileName)
            for i in xrange(max(1,-NP)):
                tl = Internal.copyRef(t2)
                bases = Internal.getNodesFromType1(tl, 'CGNSBase_t')
                for b in bases:
                    zones = Internal.getNodesFromType1(b, 'Zone_t')
                    for z in zones:
                        proc = Internal.getNodeFromName2(z, 'proc')
                        if proc is not None:
                            proc = Internal.getValue(proc)
                            if proc != i: Internal._rmNode(b, z)
                C.convertPyTree2File(tl, '%s/%s_%d.cgns'%(fileName,baseName,i))

            # Write graph
            import cPickle as pickle
            file = open('%s/graph.pk'%fileName, 'wb')
            pickle.dump(objet, file, protocol=pickle.HIGHEST_PROTOCOL)
            file.close()

#============================================================================
# Retourne le max proc pour les zones
def getMaxProc(t):
    maxProc = 0
    zones = Internal.getZones(t)
    for z in zones:
        proc = Internal.getNodeFromName2(z, 'proc')
        if proc is not None:
            proc = Internal.getValue(proc)
            maxProc = max(maxProc, proc)
    return maxProc

#==============================================================================
# Load one file
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# if graph=True,
# return the communication graph for Chimera and abutting transfers
# and the communication graph for IBM transfers
#==============================================================================
def loadFile(fileName='t.cgns', split='single', graph=False, mpirun=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if mpirun: # mpi run
        import Converter.Mpi as Cmpi
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            FILE = fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            mp = getMaxProc(t)
            if mp+1 != size: 
                raise ValueError, 'The number of mpi proc (%d) doesn t match the tree distribution (%d)'%(mp+1,size) 
            if graph:
                graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                graphIBCD = Cmpi.computeGraph(t, type='IBCD', reduction=False)
                procDict  = D2.getProcDict(t)
                procList  = D2.getProcList(t,  sort=True)
                graphN = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
            t = Cmpi.readZones(t, FILE, rank=rank)
            t = Cmpi.convert2PartialTree(t, rank=rank)
            
        else: # load 1 fichier par proc
            if graph:
                # Try to load graph from file
                if os.access('%s/graph.pk'%fileName, os.R_OK):
                    import cPickle as pickle
                    file = open('%s/graph.pk'%fileName, 'rb')
                    graphN = pickle.load(file)
                    file.close()
                # Load all skeleton proc files
                else:
                    ret = 1; no = 0; t = []
                    while ret == 1:
                       FILE = '%s_%d.cgns'%(fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          t.append(Cmpi.convertFile2SkeletonTree(FILE))
                       no += 1
                    if t != []:
                        t        = Internal.merge(t)
                        graphID  = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                        graphIBCD= Cmpi.computeGraph(t, type='IBCD', reduction=False)
                        procDict = D2.getProcDict(t)
                        procList = D2.getProcList(t,  sort=True)
                        graphN   = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                    else: print 'graph non calculable: manque de fichiers connectivite'

            FILE = '%s/%s_%d.cgns'%(fileName, baseName, rank)
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
        else: # multiple
            ret = 1; no = 0; t = []
            while ret == 1:
                FILE = '%s/%s_%d.cgns'%(fileName, baseName, no)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: t.append(C.convertFile2PyTree(FILE))
                no += 1
            if t != []: t = Internal.merge(t)
            else: t = None

    if graph: return t, graphN
    else:     return t

#==============================================================================
# Save tree in one file.
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
#==============================================================================
def saveFile(t, fileName='restart.cgns', split='single', graph=False, NP=0, mpirun=False):
    """Save tree in file."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension

    # Rip/add some useless/usefull data (FastS)
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')
    zones = Internal.getZones(t2)
    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#define')
        if node is not None:
            node = Internal.getNodeFromName1(node, 'temporal_scheme')
            if node is not None:
                integ = Internal.getValue(node)
                if integ == 'explicit':
                    C._rmVars(z, 'centers:Density_M1')
                    C._rmVars(z, 'centers:VelocityX_M1')
                    C._rmVars(z, 'centers:VelocityY_M1')
                    C._rmVars(z, 'centers:VelocityZ_M1' )
                    C._rmVars(z, 'centers:Temperature_M1')
                    C._rmVars(z, 'centers:TurbulentSANuTilde_M1')
    dtloc = Internal.getNodeFromName3(t2, '.Solver#dtloc')
    if dtloc is not None:
        node =  Internal.getNodeFromName1(t, 'TimeLevelMotion')
        if node is not None: node[1][0]= dtloc[3]
        else: Internal.createUniqueChild(t, 'TimeLevelMotion', 'DataArray_t', value=dtloc[3])
        node =  Internal.getNodeFromName1(t, 'TimeLevelTarget')
        if node is not None: node[1][0]= dtloc[4]
        else: Internal.createUniqueChild(t, 'TimeLevelTarget', 'DataArray_t', value=dtloc[4])

    # save
    if mpirun: # mpi run
        import Converter.Mpi as Cmpi
        if split == 'single':  # output in a single file
            FILE = fileName+'.cgns'
            Cmpi.convertPyTree2File(t2, FILE)
        else:
            rank = Cmpi.rank
            if rank == 0: 
                if not os.path.exists(fileName): os.makedirs(fileName)
            Cmpi.barrier()
            FILE = '%s/%s_%d.cgns'%(fileName, baseName, rank)
            C.convertPyTree2File(t2, FILE)

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            C.convertPyTree2File(t2, FILE)
        else:
            if NP == 0: NP = -getMaxProc(t2)-1
            if not os.path.exists(fileName): os.makedirs(fileName)
            for i in xrange(max(1,-NP)):
                tl = Internal.copyRef(t2)
                bases = Internal.getNodesFromType1(tl, 'CGNSBase_t')
                for b in bases:
                    zones = Internal.getNodesFromType1(b, 'Zone_t')
                    for z in zones:
                        proc = Internal.getNodeFromName2(z, 'proc')
                        if proc is not None:
                            proc = Internal.getValue(proc)
                            if proc != i: Internal._rmNode(b, z)
                C.convertPyTree2File(tl, '%s/%s_%d.cgns'%(fileName,baseName,i))

#==============================================================================
# Load t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
# dir is the directory containing files to be read 
#==============================================================================
def loadTree(fileName='t.cgns', split='single', directory='.', graph=False, mpirun=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension
    import Converter.PyTree as C

    graphN = {'graphID':None, 'graphIBCD':None, 'procDict':None, 'procList':None}
    if mpirun == True: # mpi run
        import Converter.Mpi as Cmpi
        import Distributor2.PyTree as D2
        rank = Cmpi.rank; size = Cmpi.size
        
        if split == 'single':
            # Load connect (tc)
            FILE = directory+'/'+fileName+'.cgns'
            t = Cmpi.convertFile2SkeletonTree(FILE)
            if graph:
                graphID   = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                graphIBCD = Cmpi.computeGraph(t, type='IBCD', reduction=False)
                procDict  = D2.getProcDict(t)
                procList  = D2.getProcList(t,  sort=True)
                graphN = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
            t = Cmpi.readZones(t, FILE, rank=rank)
            zones=Internal.getZones(t)
            t = Cmpi.convert2PartialTree(t, rank=rank)
            zones=Internal.getZones(t)

        else: # load 1 fichier par proc
            if graph:
                # Try to load graph
                if os.access('%s/graph.pk'%directory, os.R_OK):
                    import cPickle as pickle
                    file = open('%s/graph.pk'%directory, 'rb')
                    graphN= pickle.load(file)
                    file.close()
                # Load all sqeleton proc files
                else:
                    ret = 1; no = 0; t = []
                    while ret == 1:
                       #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, no)
                       FILE = '%s/%s%d.cgns'%(directory, fileName, no)
                       if not os.access(FILE, os.F_OK): ret = 0
                       if ret == 1: 
                          t.append(Cmpi.convertFile2SkeletonTree(FILE))
                       no += 1
                    if no == size and t != []:
                        t        = Internal.merge(t)
                        graphID  = Cmpi.computeGraph(t, type='ID'  , reduction=False)
                        graphIBCD= Cmpi.computeGraph(t, type='IBCD', reduction=False)
                        procDict = D2.getProcDict(t)
                        procList = D2.getProcList(t,  sort=True)
                        graphN   = {'graphID':graphID, 'graphIBCD':graphIBCD, 'procDict':procDict, 'procList':procList }
                    else: print 'graph non calculable: manque de fichiers connectivite'

            #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, rank)
            FILE = '%s/%s%d.cgns'%(directory, fileName, rank)
            print 'filename', FILE
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: tc= None

    else: # sequential run
        if split == 'single':
            FILE = fileName+'.cgns'
            if os.access(FILE, os.F_OK): t = C.convertFile2PyTree(FILE)
            else: t = None
        else: # multiple
            # Load connectivity (tc)
            ret = 1; no = 0; tc = []
            while ret == 1:
                #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, no)
                FILE = '%s/%s%d.cgns'%(directory, fileName, no)
                if not os.access(FILE, os.F_OK): ret = 0
                if ret == 1: tc.append(C.convertFile2PyTree(FILE))
                no += 1
            #if no == NP and t != []: t = Internal.merge(t)
            if  t != []: t = Internal.merge(t)
            else: t = None

    if graph: return t, graphN
    else:     return t

#==============================================================================
# Load t.cgns (data) tc.cgns (connectivity file) ts.cgns (stat)
# si split='single': load t.cgns
# si split='multiple': read t.cgns ou t/t_1.cgns
# autorestart possible
# return a partial tree t, a donor partial tree tc, 
# a stat partial tree,
# the communication graph for Chimera and abutting transfers
# the communication graph for IBM transfers
# dir is the directory containing files to be read 
#==============================================================================
def saveTree(t, fileName='restart.cgns', split='single', directory='.', graph=False, mpirun=False):
    """Load tree and connectivity tree."""
    import os.path
    baseName = os.path.basename(fileName)
    baseName = os.path.splitext(baseName)[0] # name without extension
    fileName = os.path.splitext(fileName)[0] # full path without extension
    import Converter.PyTree as C

    # Rip/add some useless/usefull data (FastS)
    import Converter.PyTree as C
    t2 = C.rmVars(t, 'centers:Density_P1')
    C._rmVars(t2, 'centers:VelocityX_P1')
    C._rmVars(t2, 'centers:VelocityY_P1')
    C._rmVars(t2, 'centers:VelocityZ_P1')
    C._rmVars(t2, 'centers:Temperature_P1')
    C._rmVars(t2, 'centers:TurbulentSANuTilde_P1')

    bases  = Internal.getNodesFromType1(t2    , 'CGNSBase_t')       # noeud
    own   = Internal.getNodeFromName1(bases[0], '.Solver#ownData')  # noeud
    dtloc = None
    if own is not None: dtloc = Internal.getNodeFromName1(own, '.Solver#dtloc')[1] # numpy

    zones = Internal.getZones(t2)
    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#define')
        if node is not None:
            node = Internal.getNodeFromName1(node, 'temporal_scheme')
            if node is not None:
                integ = Internal.getValue(node)
                if integ == 'explicit':
                    C._rmVars(z, 'centers:Density_M1')
                    C._rmVars(z, 'centers:VelocityX_M1')
                    C._rmVars(z, 'centers:VelocityY_M1')
                    C._rmVars(z, 'centers:VelocityZ_M1' )
                    C._rmVars(z, 'centers:Temperature_M1')
                    C._rmVars(z, 'centers:TurbulentSANuTilde_M1')

    if dtloc is not None:
       node =  Internal.getNodeFromName1(t, 'TimeLevelMotion')
       if node is not None: node[1][0]= dtloc[3]
       else: Internal.createUniqueChild(t, 'TimeLevelMotion', 'DataArray_t', value=dtloc[3])
       node =  Internal.getNodeFromName1(t, 'TimeLevelTarget')
       if node is not None: node[1][0]= dtloc[4]
       else: Internal.createUniqueChild(t, 'TimeLevelTarget', 'DataArray_t', value=dtloc[4])


    # save
    if mpirun == True: # mpi run
        import Converter.Mpi as Cmpi
        if split == 'single':  # output in a single file
            FILE = directory+'/'+fileName+'.cgns'
            Cmpi.convertPyTree2File(t2, FILE)
        else:
            rank = Cmpi.rank
            #FILE = '%s/%s_proc%d.cgns'%(directory, fileName, rank)
            FILE = '%s/%s%d.cgns'%(directory, fileName, rank)
            C.convertPyTree2File(t2, FILE)

    else: # sequential run
        FILE = directory+'/'+fileName+'.cgns'
        C.convertPyTree2File(t2, FILE)


