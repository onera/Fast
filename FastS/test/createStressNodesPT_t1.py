# - createStressNodes (pyTree) - 
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Initiator.PyTree as I
import Converter.Internal as Internal
import KCore.test as test

ni = 155 ; dx = 100./(ni-1) ; dz = 1.
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.addBC2Zone(a1, 'far', 'BCFarfield', 'imin')
a1 = C.addBC2Zone(a1, 'fam', 'FamilySpecified:mywall', 'imax')
a1 = C.fillEmptyBCWith(a1, 'wall', 'BCWall', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])
C._addFamily2Base(Internal.getNodeFromName1(t, 'Base'), 'mywall', bndType='BCWall')

# Numerics
numb = {}; numz = {}
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metric
(t, tc, metrics) = FastS.warmup(t, None)

teff = FastS.createStressNodes(t, BC=['BCWall','BCWallViscous'])
FastS._computeStress(t, teff, metrics)
Internal._rmNodesFromName(teff, 'Parameter_int')
Internal._rmNodesFromName(teff, 'Parameter_real')
test.testT(teff, 1)

teff = FastS.createStressNodes(t, BC=['mywall'])
FastS._computeStress(t, teff, metrics)
Internal._rmNodesFromName(teff, 'Parameter_int')
Internal._rmNodesFromName(teff, 'Parameter_real')
test.testT(teff, 2)
