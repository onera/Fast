# - compute (pyTree) -
# - Monobloc laminar flat plate -
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test

t = C.convertFile2PyTree('Pplane_136_96.cgns')
t = C.addState(t, MInf=0.2, ReInf=25.e6, MutSMuInf=15)

numb = { 'temporal_scheme':'implicit' }
numz = { 'scheme':'ausmpred' }
Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

FastS._applyBC(t, metrics)

teff = FastS.createStressNodes(t, BC=['BCWall'])

effort = FastS._computeStress(t, teff, metrics)

Internal._rmNodesByName(teff, '.Solver#Param')
Internal._rmNodesByName(teff, '.Solver#ownData')
test.testT(teff, 1)
