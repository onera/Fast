# - applyBC (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Initiator.PyTree as I
import KCore.Adim

MInf = 2.9
adim = KCore.Adim.adim1(MInf)

a = G.cart((0,0,0), (1,1,1), (20,20,2))
a = I.initConst(a, MInf=0.6, loc='centers')
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
a = C.addBC2Zone(a, 'extrap', 'BCExtrapolate', 'imax')
#a = C.addBC2Zone(a, 'wall', 'BCWallInviscid', 'jmin')
#a = C.addBC2Zone(a, 'farf', 'BCFarfield', 'imax')
#a = C.addBC2Zone(a, 'farf', 'BCInflowSupersonic', 'imin', data=adim)

a = C.addState(a, 'GoverningEquations', 'Euler')
a = C.addState(a, adim='adim1', MInf=MInf)
t = C.newPyTree(['Base',a])
(t, tc, metrics)  = FastS.warmup(t, None)
# Numerics
numb = {}
numb["temporal_scheme"] = "explicit"
Fast._setNum2Base(t, numb)
FastS._applyBC(t, metrics)
C.convertPyTree2File(t, 'out.cgns')
