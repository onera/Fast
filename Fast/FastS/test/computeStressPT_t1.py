# - computeStress (pyTree) - 
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Initiator.PyTree as I
import KCore.test as test

ni = 155 ; dx = 100./(ni-1) ; dz = 1.
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.fillEmptyBCWith(a1, 'wall', 'BCWall', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])

# Numerics
numb = {}; numz = {}
numb = { 'temporal_scheme': 'implicit' }
numz = { 'scheme': 'ausmpred' }
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metric
(t, tc, metrics) = FastS.warmup(t, None)

teff = FastS.createStressNodes(t, ['BCWall'])

# Compute
for nitrun in range(1,200):
    FastS._compute(t, metrics, nitrun)

effort = FastS._computeStress(t, teff, metrics)

test.testT(teff, 1)
