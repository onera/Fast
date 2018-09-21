# test en SA avec une topologie curviligne dans le plan, cartesienne en k
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Initiator.PyTree as I
import KCore.Adim as Adim
import Converter.Internal as Internal
import Dist2Walls.PyTree as DTW
import Connector.PyTree as X
import Transform.PyTree as T
import KCore.test as test

MInf = 0.1
NI = 51; NJ = 36; NK = 15
a = G.cylinder((0,0,0), 0.5, 25., 360., 0., 0.5, (NI,NJ,NK))
distrib = G.cart((0,0,0),(1/(NJ-1.),1,1),(NJ,1,1))
distrib = G.enforcePlusX(distrib,2.e-5,15,25)
a = G.map(a, distrib, dir=2)

a1 = T.subzone(a,(1,1,1),(NI/2,-1,-1)); a1[0] = 'cyl1'
a2 = T.subzone(a,(NI/2,1,1),(NI,-1,-1)); a2[0] = 'cyl2'
t = C.newPyTree(["Base"]); t[2][1][2] = [a1,a2]
t = X.connectMatch(t)

t = Internal.addGhostCells(t,t,2,2)
t = C.rmBCOfType(t, 'BCMatch')
t = C.addBC2Zone(t,'overlap','BCOverlap','imin')
t = C.addBC2Zone(t,'overlap','BCOverlap','imax')
t = C.addBC2Zone(t,'wall','BCWall','jmin')
t = C.fillEmptyBCWith(t, 'far', 'BCFarfield', dim=3)

adim = Adim.adim1(MInf=MInf, ReInf=40000., MutSMuInf=15)
t = X.applyBCOverlaps(t, depth=2)
tc = C.node2Center(t)
tc = X.setInterpData(t,tc,nature=1,loc='centers',storage='inverse',sameName=1)
tc = Internal.rmNodesByName(tc, 'GridCoordinates')
tc = Internal.rmNodesByName(tc, 'FlowSolution')
C._addState(t, MInf=MInf, alphaZ=0.)
C._addState(t, 'GoverningEquations', 'NSTurbulent')
I._initConst(t, MInf=MInf, loc='centers')
C._randomizeVar(t, 'centers:Density', 0.1, 0.1)
C._randomizeVar(t, 'centers:MomentumX', 0.01, 0.01)

# distance a la paroi
tb = C.newPyTree(['Body']); tb[2][1][2].append(G.cylinder((0,0,0), 0.5, 25., 360., 0., 0.5, (NI,1,2)))
t = DTW.distance2Walls(t, tb, loc='centers', type='ortho')

# Numerics
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 10
numb["modulo_verif"]       = 20
numz = {}
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

surf = G.cylinder((1,1,0), 0.35, 25., 360., 0., 0.5, (NI,1,2))
tsurf = C.newPyTree(['Surf']); tsurf[2][1][2]=[surf]

t    = FastS.setIBCData_zero(t, tsurf, dim=2)
#C.convertPyTree2File(t,"toto.cgns")
#import sys; sys.exit()
# Prim vars, solver tag, compact, metrics
(t,  tc, metrics) = FastS.warmup(t, tc )

# Compute
for it in xrange(1,100):
#for it in xrange(1,5):
  print 'it=',it 
  FastS._compute(t, metrics, it, tc)
  if it%20==0:
       FastS.display_temporal_criteria(t, metrics, it)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
#C.convertPyTree2File(t,'ibm.cgns')
test.testT(t,1)
