# - compute (pyTree) -
# - Monobloc laminar flat plate -
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Connector.PyTree as X
import Converter.Internal as Internal
import KCore.test as test
import sys

t = C.convertFile2PyTree('Pplane_136_96.cgns')
t = C.addState(t, MInf=0.2, ReInf=25.e6, MutSMuInf=15)

numb = { 'temporal_scheme':'implicit', 'ss_iteration':5, 'modulo_verif':50 }
numz = {'ssdom_IJK': [600, 30, 20000], 
        'io_thread': -4, 
        'time_step':0.0000152280, 
        'time_step_nature': 'local', 
        'cfl': 10., 
        'scheme':'ausmpred', 
        'epsi_newton': 0.5}
Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

FastS._applyBC(t, metrics)

nit = 500; time = 0.
for it in xrange(nit):
    FastS._compute(t, metrics, it)
    if (it%50 == 0):
        print '- %d - %g'%(it, time)
        FastS.display_temporal_criteria(t, metrics, it, format='double')
    time += numz['time_step']

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
