# - compute (pyTree) -
# - Lamb vortex [Euler/implicit] -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import KCore.test as test
import Converter.Internal as Internal

mach = 0.7
a = G.cart((0,0,0), (0.5,0.5,0.25), (200,100,2))
I._initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
t = C.newPyTree(['Base']); t[2][1][2] += [a]
C._addState(t, 'GoverningEquations', 'Euler')
C._addState(t, MInf=mach)
# Numerics
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 10
numz = {}
numz["time_step"]          = 0.03
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

nit = 1000; time = 0.
timeStep = numz['time_step']
for it in xrange(nit):
    FastS._compute(t, metrics, it)
    time += timeStep

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
