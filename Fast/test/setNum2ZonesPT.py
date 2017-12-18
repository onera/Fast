# - setNum2Zones (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Fast.PyTree as Fast
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base', a])

numz = { 'scheme': 'ausmpred',
         'time_step_nature': 'global',
         'time_step': 0.001 }

Fast._setNum2Zones(t, numz)
Internal.printTree(t)
