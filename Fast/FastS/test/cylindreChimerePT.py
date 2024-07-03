# - FastS: Cylindre Chimere -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import numpy
import Transform.PyTree as T
import Initiator.PyTree as I
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Converter.Internal as Internal
import KCore.test as test
import sys

MInf = 0.5; alpha = 0.

# Solver settings
numb = {'temporal_scheme':'implicit', 'ss_iteration':5, 'modulo_verif':7} 
numz = {'time_step':0.01, 'scheme':'senseur', 'extract_res':2}

# Grille cartesienne
NI = 501; dh = 40./(NI-1)
a = G.cart((-20.,-20.,0),(dh,dh,1.),(NI,NI,2))
a = C.fillEmptyBCWith(a, 'far', 'BCFarfield', dim=2) 

# Cylindre
b = G.cylinder((0,0,0), 1., 2.0, 180., 0., 5*0.01, (80,20,1))
c = G.cylinder((0,0,0), 1., 2.0, 360., 180., 5*0.01, (80,20,1))
t = C.newPyTree(['Base']); t[2][1][2] += [b,c]
t = X.connectMatch(t, dim=2)
t = Internal.addGhostCells(t, t, 2, adaptBCs=1)
t = C.rmBCOfType(t, 'BCMatch')
t = T.addkplane(t)

zones = Internal.getNodesFromType(t, 'Zone_t')
b = zones[0]; c = zones[1]
b = C.addBC2Zone(b, 'wall', 'BCWall', 'jmin')
b = C.addBC2Zone(b, 'ov', 'BCOverlap', 'jmax')
b = C.addBC2Zone(b, 'ov', 'BCOverlap', 'imin')
b = C.addBC2Zone(b, 'ov', 'BCOverlap', 'imax')

c = C.addBC2Zone(c, 'wall', 'BCWall', 'jmin')
c = C.addBC2Zone(c, 'ov', 'BCOverlap', 'jmax')
c = C.addBC2Zone(c, 'ov', 'BCOverlap', 'imin')
c = C.addBC2Zone(c, 'ov', 'BCOverlap', 'imax')

t = C.newPyTree(['Cylindre', 'Cart']) 
t[2][1][2] += [b,c]; t[2][2][2] += [a]
t = C.fillMissingVariables(t)
#C.convertPyTree2File(t, 'out.cgns'); sys.exit()

# Body
s = D.circle((0,0,0), 1.)
s = T.addkplane(s, 2); s = T.translate(s, (0,0,-0.1))
s = C.convertArray2Tetra(s); s = G.close(s)
e = P.exteriorFaces(s); e = T.splitConnexity(e)
p = G.fittingPlaster(e[0], bumpFactor=0.) 
g1 = G.gapfixer(e[0], p)
p = G.fittingPlaster(e[1], bumpFactor=0.) 
g2 = G.gapfixer(e[1], p)
s = T.join([s,g1,g2]); s = G.close(s)

# Blanking
bodies = [[s]]
BM = numpy.array([[0],[1]], Internal.E_NpyInt)
t = X.applyBCOverlaps(t, depth=2)
t = X.blankCells(t, bodies, BM, blankingType='cell_intersect')
t = X.setHoleInterpolatedPoints(t, depth=+2)

tc = C.node2Center(t)
tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
                     sameName=1, method='leastsquares', dim=2)
tc = C.rmVars(tc, 'cellN') # tres important pour l'instant

# Init
t = C.addState(t, 'GoverningEquations', 'NSLaminar')
t = C.addState(t, MInf=MInf, alphaZ=alpha)
t = I.initConst(t, MInf=MInf, alphaZ=alpha, loc='centers')
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, tc)

#sys.exit()
nit = 585
#nit = 1

for it in range(nit):
    print('it=' ,it)
    FastS._compute(t, metrics, it, tc)
    if it%14==0:
       FastS.displayTemporalCriteria(t, metrics, it)

C.convertPyTree2File(t,'verif.cgns')

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
#test.testT(t, 1)
