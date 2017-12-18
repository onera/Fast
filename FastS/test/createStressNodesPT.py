# - computeStress (pyTree) - 
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Initiator.PyTree as I

ni = 155 ; dx = 100./(ni-1) ; dz = 1.
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.fillEmptyBCWith(a1, 'wall', 'BCWall', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])

# Numerics
numb = {}; numz = {}
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metric
(t, tc, metrics) = FastS.warmup(t, None)

teff = FastS.createStressNodes(t, BC=['BCWall','BCWallViscous'])

C.convertPyTree2File(teff, 'stress.cgns')
