# - compute (pyTree) -
# - Lamb vortex on octree -
import Connector.PyTree as X
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import FastS.Mpi as FastSmpi
import Converter.PyTree as C
import Distributor2.PyTree as D2
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import sys

MInf = 0.7
rank = Cmpi.rank ; size = Cmpi.size
FILE = 'lamb.cgns' ; FILED = 'lambD.cgns'

# lecture du squelette
t = Cmpi.convertFile2SkeletonTree(FILE)
tc = Cmpi.convertFile2SkeletonTree(FILED)

# equilibrage
(t, dic) = D2.distribute(t, NProc=size, algorithm='fast', useCom=0)
tc = D2.copyDistribution(tc, t)
graph = Cmpi.computeGraph(tc, type='ID')
procDict = D2.getProcDict(tc)

# load des zones locales dans le squelette
t = Cmpi.readZones(t, FILE, rank=rank)
tc = Cmpi.readZones(tc, FILED, rank=rank)

t = Cmpi.convert2PartialTree(t)
tc = Cmpi.convert2PartialTree(tc)
# Init
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, MInf=MInf)
# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numz = {}
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, tc)

graph = {'graphID':graph, 'graphIBCD':None, 'procDict':procDict,  'procList':[0] }

nit = 1000 ; time = 0.
for it in range(nit):
    FastSmpi._compute(t, metrics, it, tc, graph)
    if rank == 0 and it%10 == 0:
        print('- %d - %g -'%(it, time)); sys.stdout.flush()
    time += numz['time_step']
Cmpi.convertPyTree2File(t, 'out.cgns')
