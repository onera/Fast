# - compute (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import Transform.PyTree as T
import Fast.PyTree as Fast

num = {'dim':2, 'model':'Euler',
       'integ':'unsteady', 'timeStep':0.01, 'CFL':0.1,
       'scheme':'ausmup'}
numC = num.copy(); numC['solver'] = 'FastS'
numIJK = num.copy(); numIJK['solver'] = 'FastS'
numP = num.copy(); numP['solver'] = 'FastP'

# Grille cartesienne (FastS)
a = G.cart((-3.5,0,0), (0.05,0.05,0.25), (50,50,6))
a = C.addBC2Zone(a, 'left', 'BCFarfield', 'imin')
a = C.addBC2Zone(a, 'right', 'BCExtrapolate', 'imax')
a = C.addBC2Zone(a, 'wall', 'BCWallInviscid', 'jmin')
a = C.addBC2Zone(a, 'up', 'BCExtrapolate', 'jmax')
Fast._setNum2Zones(a, numC)

# Grille curviligne (FastS)
b = G.cylinder((0.2,0,0), 0.5, 2.5, 150., 0., 1.25, (30,30,6))
b = C.addBC2Zone(b, 'left', 'BCFarfield', 'imin')
b = C.addBC2Zone(b, 'right', 'BCExtrapolate', 'imax')
b = C.addBC2Zone(b, 'wall', 'BCExtrapolate', 'jmin')
b = C.addBC2Zone(b, 'up', 'BCExtrapolate', 'jmax')
Fast._setNum2Zones(b, numIJK)

# Grille polyedrique (FastP)
c = G.cartTetra((2.2,0,0), (0.05,0.05,0.25), (50,50,6))
c = T.dual(c)
Fast._setNum2Zones(c, numP)

t = C.newPyTree(['Base',a,b,c])
t = C.addState(t, MInf=0.5, alphaZ=-1., alphaY=0.)
t = I.initConst(t, MInf=0.5, alphaZ=-1., alphaY=0., loc='centers')

t, metrics = Fast.metric(t)
t = C.fillMissingVariables(t)
#C.convertPyTree2File(t, 'out.cgns'); import sys; sys.exit()

nit = 100; time = 0.
timeStep = Fast.getValueFromTag(t, 'timeStep')
for it in xrange(nit):
    t = Fast.applyBC(t, metrics)
    t = Fast.compute(t, metrics)
    if it%10 == 0:
        print '- %d - %g'%(it, time)
        CPlot.display(t, dim=2, mode=3)
    time += timeStep

t = Fast.applyBC(t, metrics)
C.convertPyTree2File(t, 'out.cgns')
