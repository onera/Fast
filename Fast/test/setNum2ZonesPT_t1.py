# - setNum2Zones (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Generator.PyTree as G
import KCore.test as test

t = C.newPyTree(['Base'])
a = G.cart((0,0,0), (1,1,1), (10,10,10))
t[2][1][2] += [a]

numz = {'scheme': 'ausmpred',
        'time_step_nature': 'global',
        'time_step': 0.001}

Fast._setNum2Zones(t, numz)
test.testT(t, 1)
