# - compute (pyTree) -
# - reflexion d'un choc sur une plaque plane -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Initiator.PyTree as I
import Fast.PyTree as Fast
import Transform.PyTree as T
import KCore.Adim as Adim
import KCore.test as test

MInf = 2.9
adim = Adim.adim1(MInf)
state1 = [adim[0], adim[1], adim[2], adim[3], adim[4]]
ro2 = (1.52819/adim[5])**(1./1.4)*adim[0]
u2 = adim[1]/adim[0]; v2 = adim[2]/adim[0]; w2 = adim[3]/adim[0]
state2 = [1.6999600, 4.4527732, -0.86072205, 0., 9.8700385]

N = 101; h = 4.1/(N-1)
a = G.cart((0,0,0), (h,h,0.001), (N,int(N/4.1),1))
a = Internal.addGhostCells(a, a, 2, adaptBCs=0)
a = T.addkplane(a)
C._addBC2Zone(a, 'wall', 'BCWallInviscid', 'jmin')
#C._addBC2Zone(a, 'super', 'BCFarfield', 'imin', data=state1)
#C._addBC2Zone(a, 'super', 'BCFarfield', 'jmax', data=state2)
C._addBC2Zone(a, 'super', 'BCInflowSupersonic', 'jmax', data=state2)
C._addBC2Zone(a, 'super', 'BCInflowSupersonic', 'imin', data=state1)
C._addBC2Zone(a, 'extr', 'BCExtrapolate', 'imax')
t = C.newPyTree(['Base', a])
I._initConst(t, MInf=MInf, loc='centers')
C._addState(t, 'GoverningEquations', 'Euler')
C._addState(t, MInf=MInf)

numb = {'temporal_scheme':'explicit', 'ss_iteration':20} 
numz = {'time_step':0.0007*2., 'scheme':'ausmpred'}
Fast._setNum2Base(t, numb)
Fast._setNum2Zones(t, numz);
#C.convertPyTree2File(t,'test.cgns')
#(t, tc, metrics) = FastS.warmup(t, None)
(t, tc, metrics) = Fast.warmup(t, None)

#nit = 2300
nit = 23
#nit=1
for it in range(nit):
    Fast._compute(t, metrics, it, NIT=100)
    #Fast._compute(t, metrics, it, NIT=1)


Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
