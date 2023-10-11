# - compute (pyTree) -
# - Lamb vortex sur 2 domaines -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Post.PyTree as P
import Connector.PyTree as X

a1 = G.cart((0,0,0), (0.5,0.5,0.25), (105,100,2))
a1 = C.addBC2Zone(a1, 'ov', 'BCOverlap', 'imax')
a2 = G.cart((50,0,0), (0.5,0.5,0.25), (100,100,2))
a2 = C.addBC2Zone(a2, 'ov', 'BCOverlap', 'imin')
t = C.newPyTree(['Base',[a1,a2]])

# Chimera
t = X.applyBCOverlaps(t, depth=2)
tc = C.node2Center(t)
tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
                     sameName=1, method='leastsquares', dim=2)
tc = C.rmVars(tc, 'cellN') # tres important pour l'instant

# Init
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, MInf=0.7)
t = I.initLamb(t, position=(40.,25.), Gamma=2., MInf=0.7, loc='centers')

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 4
numb["modulo_verif"]       = 5
numz = {}
numz["time_step_nature"]          = 'local'
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
FastC._setNum2Zones(t, numz) ; FastC._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

nit = 11 ; time = 0.
timeStep = numz['time_step']
for it in range(nit):
    FastS._compute(t, metrics, it, tc)
    if it%5 == 0:
        print('- %d - %g'%(it, time))
        FastS.display_temporal_criteria(t,metrics,it)
        #CPlot.display(t, dim=2, mode=3, scalarField=1)
    time += timeStep
