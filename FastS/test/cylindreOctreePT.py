# - FastS: Cylindre Chimere Octree -

import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import numpy
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import CPlot.PyTree as CPlot
import Converter.Internal as Internal
import KCore.Adim as Adim
import sys

MInf = 0.5 ; alpha = 0.
adim = Adim.adim1(MInf=MInf)

# Cylindre
b = G.cylinder((0,0,0), 1., 2.0, 180., 0., 5*0.01, (80,20,1))
c = G.cylinder((0,0,0), 1., 2.0, 360., 180., 5*0.01, (80,20,1))
t = C.newPyTree(['Base']) ; t[2][1][2] += [b,c]
t = X.connectMatch(t, dim=2)
t = Internal.addGhostCells(t, t, 2, adaptBCs=1)
t = C.rmBCOfType(t, 'BCMatch')
t = T.addkplane(t)

zones = Internal.getZones(t)
b = zones[0]; c = zones[1]
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'ov', 'BCOverlap', 'jmax')
b = C.addBC2Zone(b, 'ov', 'BCOverlap', 'imin')
b = C.addBC2Zone(b, 'ov', 'BCOverlap', 'imax')

c = C.addBC2Zone(c, 'wall', 'BCWall', 'jmin')
c = C.addBC2Zone(c, 'ov', 'BCOverlap', 'jmax')
c = C.addBC2Zone(c, 'ov', 'BCOverlap', 'imin')
c = C.addBC2Zone(c, 'ov', 'BCOverlap', 'imax')

# Octree avec grilles cartesiennes
NI = 501 ; dh = 40./(NI-1)
vmin = 31
ovs1 = C.extractBCOfType([b],"BCOverlap")[0] 
dims = Internal.getZoneDim(ovs1); ni = dims[1]; nj = 1
ovs1 = T.subzone(ovs1,(1,1,1),(ni,nj,1))

ovs2 = C.extractBCOfType([c],"BCOverlap")[0] 
dims = Internal.getZoneDim(ovs2); ni = dims[1]; nj = 1
ovs2 = T.subzone(ovs2,(1,1,1),(ni,nj,1))

snears = [dh*(vmin-1)]*2
o = G.octree([ovs1,ovs2], snears, dfar=100., balancing=1)

res = G.octree2Struct(o,vmin=vmin,ext=3,optimized=0,merged=1)
res = C.fillEmptyBCWith(res, 'far', 'BCFarfield', dim=2) 
res = T.addkplane(res)

t = C.newPyTree(['Cylindre', 'Cart']) 
t[2][1][2] += [b,c] ; t[2][2][2] += res
t = C.fillMissingVariables(t)

# Body
s = D.circle((0,0,0), 1.)
s = T.addkplane(s, 2) ; s = T.translate(s, (0,0,-0.1))
s = C.convertArray2Tetra(s) ; s = G.close(s)
e = P.exteriorFaces(s) ; e = T.splitConnexity(e)
p = G.fittingPlaster(e[0], bumpFactor=0.) 
g1 = G.gapfixer(e[0], p)
p = G.fittingPlaster(e[1], bumpFactor=0.) 
g2 = G.gapfixer(e[1], p)
s = T.join([s,g1,g2]); s = G.close(s)

# Blanking
bodies = [[s]]
BM = numpy.array([[0],[1]],numpy.int32)
t = X.applyBCOverlaps(t, depth=2)
t = X.blankCells(t, bodies, BM, blankingType='cell_intersect')
t = X.setHoleInterpolatedPoints(t, depth=+2)
t = X.optimizeOverlap(t,priorities=['Cylindre',0,'Cart',1])
t = X.maximizeBlankedCells(t,depth=2)
tc = C.node2Center(t)
tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
                     sameName=1, method='leastsquares',dim=2)
# Init
t = C.addState(t, MInf=MInf, alphaZ=alpha)
t = C.addState(t, 'GoverningEquations', 'Euler')
# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numz = {}
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)


t = I.initConst(t, MInf=MInf, alphaZ=alpha, loc='centers')
(t, tc, metrics)  = FastS.warmup(t, tc)

nit = 1000
for it in range(nit):
    FastS._compute(t, metrics, nit, tc)
    if it%100 == 0:
        print('- %d -'%it); sys.stdout.flush()
        CPlot.display(t, dim=2, mode=3, scalarField=1)
C.convertPyTree2File(t, "out.cgns")
