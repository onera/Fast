# - compute (pyTree) -
# - Lamb vortex -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import KCore.test as test
import Converter.Internal as Internal

mach = 0.7
##
## Dim test cpu (28.7 versus 30.6s pour 300dt 6thread)
##
a = G.cart((0,0,0), (0.5,0.5,0.25), (201,101,32))
I._initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
t = C.newPyTree(['Base', a])

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numz = {}
#numz["omp"]                = "scater"
#numz["io_thread"]          = 11
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
t = C.addState(t, 'GoverningEquations', 'NSLaminar')
t = C.addState(t, MInf=mach)
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

nit = 50; time = 0.
timeStep = numz['time_step']
for it in range(nit):
    FastS._compute(t, metrics, it)
    if it%100 == 0:
        print('- %d - %g'%(it, time))
    time += timeStep

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
