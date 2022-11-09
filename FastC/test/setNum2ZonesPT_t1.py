# - setNum2Zones (pyTree) -
import Converter.PyTree as C
import FastC.PyTree as FastC
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base', a])

numz = {'scheme': 'ausmpred',
        'time_step_nature': 'global',
        'time_step': 0.001}

FastC._setNum2Zones(t, numz)
test.testT(t, 1)
