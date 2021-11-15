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
test.TOLERANCE = 1.e-10

MInf = 0.1
NI = 51; NJ = 36; NK = 15
a = G.cylinder((0,0,0), 0.5, 25., 360., 0., 0.5, (NI,NJ,NK))
distrib = G.cart((0,0,0),(1/(NJ-1.),1,1),(NJ,1,1))
distrib = G.enforcePlusX(distrib,2.e-5,15,25)
a = G.map(a,distrib,dir=2)

a1 = T.subzone(a,(1,1,1),(NI//2,-1,-1)); a1[0] = 'cyl1'
a2 = T.subzone(a,(NI//2,1,1),(NI,-1,-1)); a2[0] = 'cyl2'
t = C.newPyTree(["Base",a1,a2])
t = X.connectMatch(t)
Internal._addGhostCells(t,t,2,adaptBCs=0)
C._rmBCOfType(t, 'BCMatch')
C._addBC2Zone(t,'overlap','BCOverlap','imin')
C._addBC2Zone(t,'overlap','BCOverlap','imax')
C._addBC2Zone(t,'wall','BCWall','jmin')
C._fillEmptyBCWith(t, 'far', 'BCFarfield', dim=3)

adim = Adim.adim1(MInf=MInf, ReInf=40000., MutSMuInf=15)
X._applyBCOverlaps(t, depth=2)
tc = C.node2Center(t)
tc = X.setInterpData(t,tc,nature=1,loc='centers',storage='inverse',sameName=1)
Internal._rmNodesByName(tc,'GridCoordinates')
Internal._rmNodesByName(tc,'FlowSolution')
C._addState(t, MInf=MInf, alphaZ=0.)
C._addState(t, 'GoverningEquations', 'NSTurbulent')
I._initConst(t, MInf=MInf, loc='centers')
C._randomizeVar(t,'centers:Density',0.1,0.1)
C._randomizeVar(t,'centers:MomentumX',0.01,0.01)

# distance a la paroi
tb = C.newPyTree(['Body']); tb[2][1][2].append(G.cylinder((0,0,0), 0.5, 25., 360., 0., 0.5, (NI,1,2)))
DTW._distance2Walls(t,tb,loc='centers',type='ortho')
#
# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numz = {}
numz["time_step"]          = 0.0001
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metrics
(t, tc, metrics) = FastS.warmup(t, tc)

#FastS._applyBC(t, metrics)

teff = FastS.createStressNodes(t, BC=['BCWall'])

for it in range(1,200): FastS._compute(t, metrics, it, tc)

effort = FastS._computeStress(t, teff, metrics)
#on supprime ghost de l'arbre effort car valeur non initialisee
Internal._rmGhostCells(teff,teff,2,2)
Internal._rmNodesByName(teff, '.Solver#Param')
Internal._rmNodesByName(teff, '.Solver#ownData')
test.testT(teff, 1)
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 2)
