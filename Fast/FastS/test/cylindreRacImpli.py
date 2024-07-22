# test en SA avec une topologie curviligne dans le plan, cartesienne en k
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Initiator.PyTree as I
import KCore.Adim as Adim
import Converter.Internal as Internal
import Dist2Walls.PyTree as DTW
import Connector.PyTree as X
import Transform.PyTree as T

MInf = 0.1
NI = 51; NJ = 36; NK = 15
a = G.cylinder((0,0,0), 0.5, 25., 360., 0., 0.5, (NI,NJ,NK))
distrib = G.cart((0,0,0),(1/(NJ-1.),1,1),(NJ,1,1))
distrib = G.enforcePlusX(distrib,2.e-5,15,25)
a = G.map(a, distrib, dir=2)

a1 = T.subzone(a,(1,1,1),(int(NI/2),-1,-1)); a1[0] = 'cyl1'
a2 = T.subzone(a,(int(NI/2),1,1),(NI,-1,-1)); a2[0] = 'cyl2'
t = C.newPyTree(["Base",a1,a2])
t = X.connectMatch(t)

#duplication raccord match en BC classique  (pour implicitation)
zones = Internal.getZones(t)
for z in zones:
  connect =  Internal.getNodeFromName3(z, 'ZoneGridConnectivity')
  matchs  =  Internal.getNodesFromType1(connect, 'GridConnectivity1to1_t')
  for match in matchs:
     rg =  Internal.getNodeFromName1(match, 'PointRange')[1]
     win=[rg[0,0],rg[0,1],rg[1,0],rg[1,1], rg[2,0],rg[2,1] ] 
     C._addBC2Zone(z, 'racinf', 'BCRacinf', win)

C._addBC2Zone(t,'wall','BCWall','jmin')
C._fillEmptyBCWith(t, 'far', 'BCFarfield', dim=3)


Internal._addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
C._rmBCOfType(t, 'BCMatch')
C._addBC2Zone(t,'overlap','BCOverlap','imin')
C._addBC2Zone(t,'overlap','BCOverlap','imax')
C._fillEmptyBCWith(t, 'far', 'BCFarfield', dim=3)

adim = Adim.adim1(MInf=MInf, ReInf=40000., MutSMuInf=15)
X._applyBCOverlaps(t, depth=2)
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
DTW._distance2Walls(t, tb, loc='centers', type='ortho')

C._rmVars(t, 'centers:cellN')

# Numerics
modulo_verif = 50
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 10
numb["modulo_verif"]       = modulo_verif
numz = {}
numz["time_step"]          = 0.001
numz["scheme"]             = "ausmpred"
numz["lu_match"]           = 1
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

# Prim vars, solver tag, compact, metrics
(t,  tc, metrics) = FastS.warmup(t, tc )

# Compute
for it in range(1,40): 
   FastS._compute(t, metrics, it, tc)
   if it%modulo_verif == 0:
        print('it=%d'%it)
        FastS.display_temporal_criteria(t, metrics, it, format='double')

C.convertPyTree2File(t, 'out.cgns')