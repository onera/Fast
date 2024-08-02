# - computeStress for flow rate (pyTree) - 
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Converter.Internal as Internal
import Initiator.PyTree as I
import numpy
import KCore.test as test

ni = 155 ; dx = 100./(ni-1) ; dz = 1.
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.fillEmptyBCWith(a1, 'inflow', 'BCInflow', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])

# Numerics
numb = {}; numz = {}
numb = { 'temporal_scheme':'implicit' }
numz = { 'scheme':'ausmpred' }
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metric
(t, tc, metrics) = FastS.warmup(t, None)

flowRateTree = FastS.createStressNodes(t, BC=['BCInflow'])

# Compute
for nitrun in range(1,200):
    FastS._compute(t, metrics, nitrun)

effort = FastS._computeStress(t, flowRateTree, metrics)

# compute of the mass flow rate with numerical flux old fashion
flowRateZones = Internal.getZones(flowRateTree)
flowRate = 0.
for z in flowRateZones:
    sol      = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
    density  = Internal.getNodeFromName1(sol, 'Density')[1]
    flowRate  += numpy.sum(density)

print('the mass flow rate accross the Inflow BC is: ', flowRate)
test.testO(flowRate, 1)

# compute of the mass flow rate with numerical flux new way
print('the mass flow rate accross the Inflow BC is: ', effort[7])
test.testO(effort[7], 2)
