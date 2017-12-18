# - compute (pyTree) -
# - Lamb vortex on octree -
import Connector.PyTree as X
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import FastS.Mpi as FastSmpi
import Converter.PyTree as C
import Distributor2.PyTree as D2
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import KCore.test as test
import sys

#mpirun -np 2 -genv OMP_NUM_THREADS=1 python lambOctreeDPT_t1.py
MInf = 0.7

rank = Cmpi.rank; size = Cmpi.size
FILE = 'lamb.cgns'
FILED = 'lambD.cgns'
#FILED = 'tc_compact.cgns'
# lecture du squelette
t = Cmpi.convertFile2SkeletonTree(FILE)
tc = Cmpi.convertFile2SkeletonTree(FILED)

# equilibrage
print 'size mpi', size
(t, dic) = D2.distribute(t, NProc=size, algorithm='fast', useCom=0)
tc = D2.copyDistribution(tc, t)
graph = Cmpi.computeGraph(tc, type='ID')
procDict = D2.getProcDict(tc)

print 'dic', procDict


# load des zones locales dans le squelette
t = Cmpi.readZones(t, FILE, rank=rank)
tc = Cmpi.readZones(tc, FILED, rank=rank)

t = Cmpi.convert2PartialTree(t)
tc = Cmpi.convert2PartialTree(tc)

Cmpi.convertPyTree2File(t, 't1.cgns')
Cmpi.convertPyTree2File(tc, 't1c.cgns')
#sys.exit()
Cmpi.barrier()
t,tc,ts,graph=Fast.load('t1.cgns','t1c.cgns', split='single', restart=False , NP=size)

print 'graph', graph['graphID']
# Init
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, MInf=MInf)
# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numz = {}
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

#Initialisation parametre calcul: calcul metric + var primitive + compactage + alignement + placement DRAM
graph1 ={'graphID':graph, 'graphIBCD':None, 'procDict':procDict}
(t, tc, metrics) = FastS.warmup(t, tc, graph=graph)

if rank==0:
  print 'graph', graph['graphID'],  graph['graphIBCD']
  print 'dict ', procDict
#C.convertPyTree2File(t , 't_test'+str(rank)+'.cgns')
#C.convertPyTree2File(tc, 'tc_test'+str(rank)+'.cgns')
#sys.exit()

nit = 1000; time = 0.
for it in xrange(nit):
    #FastSmpi._compute(t, metrics, it, tc, graph)
    FastS._compute(t, metrics, it, tc)
    if (rank == 0 and it%100 == 0):
        print '- %d - %g -'%(it, time); sys.stdout.flush()
    time += numz['time_step']

#Cmpi.convertPyTree2File(t, 'out.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
