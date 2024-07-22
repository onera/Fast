# - displayTemporalCriteria (pyTree)-
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Initiator.PyTree as I
import KCore.test as test

ni = 155; dx = 100./(ni-1); dz = 0.01
a = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a = C.fillEmptyBCWith(a, 'far', 'BCFarfield', dim=2)
I._initConst(a, MInf=0.4, loc='centers')
C._addState(a, 'GoverningEquations', 'Euler')
C._addState(a, MInf=0.4)
t = C.newPyTree(['Base',a])

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numz = {}
numz["time_step"]          = 0.00004444
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics)  = FastS.warmup(t, None)

# Compute
for it in range(1,10):
    FastS._compute(t, metrics, it)
    FastS._displayTemporalCriteria(t, metrics, it, format='double')

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
