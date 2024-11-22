# - compute (pyTree) -
# - Lamb vortex sur 2 domaines -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Connector.PyTree as X
import KCore.Adim as Adim
import KCore.test as test

MInf = 0.7
adim = Adim.adim1(MInf=MInf)
a1 = G.cart((0,0,0), (0.5,0.5,0.25), (105,100,2))
a1 = C.addBC2Zone(a1, 'ov', 'BCOverlap', 'imax')
a2 = G.cart((50,0,0), (0.5,0.5,0.25), (100,100,2))
a2 = C.addBC2Zone(a2, 'ov', 'BCOverlap', 'imin')
t = C.newPyTree(['Base']); t[2][1][2] += [a1,a2]

# Chimera
t = X.applyBCOverlaps(t, depth=2)
tc = C.node2Center(t)
tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse',
                     sameName=1, method='leastsquares', dim=2)

t = I.initLamb(t, position=(45.,25.), Gamma=2., MInf=0.7, loc='centers')
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, MInf=MInf)
# Numerics
numb = {'temporal_scheme':'explicit', 'ss_iteration':20}
numz = {'time_step':0.01, 'scheme':'ausmpred'}
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)
#Initialisation parametre calcul: calcul metric + var primitive + compactage + alignement + placement DRAM
(t, tc, metrics)  = FastS.warmup(t, tc)

nit = 1000; time = 0.
timeStep = numz['time_step']
for it in range(nit):
    FastS._compute(t, metrics, it, tc)
    time += timeStep

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
