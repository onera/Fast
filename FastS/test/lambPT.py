# - compute (pyTree) -
# - Lamb vortex -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Post.PyTree as P
import Converter.Internal as Internal
import sys

mach = 0.7
a = G.cart((0,0,0), (0.5,0.5,0.25), (200,100,2))
I._initLamb(a, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')
t = C.newPyTree(['Base', a])
C._addState(t, 'GoverningEquations', 'Euler')
C._addState(t, MInf=mach)

# Numerics
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 10
numb["modulo_verif"]       = 5
numz = {}
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

zones      = Internal.getNodesFromType2(t, 'Zone_t')
param_int  = Internal.getNodeFromName2(zones[0], 'Parameter_int')[1]
param_int[31] = 150

nit = 10 ; time = 0.
timeStep = numz['time_step']
for it in xrange(nit):
    FastS._compute(t, metrics, it)
    if (it%1 == 0):
        print '- %d - %g'%(it, time)
        #CPlot.display(t, dim=2, mode=3, isoEdges=1)
    time += timeStep
