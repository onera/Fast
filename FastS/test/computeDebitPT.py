# - computeStress (pyTree) - 
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Converter.Internal as Internal
import Initiator.PyTree as I
import numpy

ni = 155 ; dx = 100./(ni-1) ; dz = 1.
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.fillEmptyBCWith(a1, 'wall', 'BCInflow', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])

# Numerics
numb = {}; numz = {}
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metric
(t, tc, metrics) = FastS.warmup(t, None)

debit_inflow = FastS.createStressNodes(t, BC=['BCInflow'])


# Compute
for nitrun in range(1,200):
    FastS._compute(t, metrics, nitrun)

effort = FastS._computeStress(t, debit_inflow , metrics)

# compute of the mass flow rate with numerical flux old fashion
zones_debit = Internal.getZones(debit_inflow)
debit =0.
for z in zones_debit:
    sol      = Internal.getNodesFromName1(z,'FlowSolution#Centers')
    density  = Internal.getNodesFromName1(sol,'Density')[1] 
    debit   +=  numpy.sum(density)

print 'the mass flow rate accross the Inflow BC is', debit

# compute of the mass flow rate with numerical flux new fashion
print 'the mass flow rate accross the Inflow BC is', zones_debit[7]

C.convertPyTree2File(debit_inflow, 'debit.cgns')
