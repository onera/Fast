# - metric (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T
import Fast.PyTree as Fast

numC = {'solver': 'FastS'}
numIJK = {'solver': 'FastS'}
numP = {'solver': 'FastP'}

# Grille cartesienne (FastS)
a = G.cart((-3.5,0,0), (0.05,0.05,0.25), (50,50,6))
Fast._setNum2Zones(a, numC)

# Grille curviligne (FastS)
b = G.cylinder((0.2,0,0), 0.5, 2.5, 150., 0., 1.25, (30,30,6))
Fast._setNum2Zones(b, numIJK)

# Grille polyedrique (FastP)
c = G.cartTetra((2.2,0,0), (0.05,0.05,0.25), (50,50,6))
c = T.dual(c)
Fast._setNum2Zones(c, numP)

t = C.newPyTree(['Base',a,b,c])

t, metrics = Fast.metric(t)
C.convertPyTree2File(t, 'out.cgns')
