# - computeVariables (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Initiator.PyTree as I
import KCore.test as test

ni = 155; dx = 100./(ni-1); dz = 0.01
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.fillEmptyBCWith(a1, 'far', 'BCFarfield', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, 'EquationDimension', 2)
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numz = {}
numz["time_step"]          = 0.00004444
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metric
(t, tc, metrics) = FastS.warmup(t, None)

C._initVars(t, 'centers:Enstrophy', 0.)

# Compute
for it in range(1,200):
    FastS._compute(t, metrics, it)
    FastS._computeVariables(t, metrics, ['Enstrophy'])

# Save
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

#on supprime ghost de l'arbre effort car valeur non initialisee
t = FastS.rmGhostCells(t,2,1)
test.testT(t, 1)
