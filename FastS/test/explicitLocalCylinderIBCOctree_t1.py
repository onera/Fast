# - FastS: Cylindre IBC Octree -
# - This is based on cylindreIBCOctreePT.py -
# - explicit local or Local time stepping (LTS) -
import CPlot.PyTree as CPlot
import Connector.Connector as Connector
import Connector.OversetData as OversetData
import Connector.PyTree as X
import Converter
import Converter.Internal as Internal
import Converter.PyTree as C
import Dist2Walls.PyTree as DTW
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Generator.PyTree as G
import Geom.PyTree as D
import Initiator.PyTree as I
import KCore
import KCore.Adim as Adim
import KCore.test as test
import Post.PyTree as P
import Transform.PyTree as T
import numpy
import sys

explicit_select = "explicit_local"

MInf = 0.2 ; alpha = 0.
adim = Adim.adim1(MInf=MInf)

# Cylindre
s = G.cylinder((0,0,0), 1., 2.0, 360., 0., 5*0.01, (160,1,1))

# Grilles cartesiennes
NI = 501 ; dh = 40./(NI-1); vmin = 31
snears = [dh*(vmin-1)]
o      = G.octree([s], snears,dfar = 10., balancing=1)
res    = G.octree2Struct(o,vmin=vmin,ext=3,optimized=0,merged=1)
res    = C.fillEmptyBCWith(res, 'far', 'BCFarfield', dim=2)
res    = T.addkplane(res)

t = C.newPyTree(['Cart']); t[2][1][2] = res
t = C.fillMissingVariables(t)

# Interpolations pour Octree
t = X.applyBCOverlaps(t, depth=2)
tc = C.node2Center(t)
if explicit_select != "explicit_local":
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse',sameName=1, method='lagrangian',dim=2)

t1c= Internal.copyRef(tc)
test.testT(t1c, 4)

# Init
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, 'EquationDimension', 2)
t = C.addState(t, MInf=MInf, alphaZ=alpha)
t = I.initConst(t, MInf=MInf, alphaZ=alpha, loc='centers')
dt = 0.01;
time_lvl_max = 0
list_lvl2=['cart3','cart4','cart6','cart8']
if explicit_select == "explicit_local":
    zones = Internal.getNodesFromType(t, 'Zone_t')
    for z in zones:
       if z[0]=="cart5":
           z          = C._initVars(z,'centers:CFL',3)
       elif z[0] in list_lvl2:
           z          = C._initVars(z,'centers:CFL',2)
       else:
           z          = C._initVars(z,'centers:CFL',1)
    for z in zones:
       dim        = Internal.getZoneDim(z)
       zp         = T.subzone(z, (3,3,1), (dim[1]-2,dim[2]-2,1))
       cflmax_loc = C.getMaxValue(zp, 'centers:CFL')
       niveauTps  = pow(2,(cflmax_loc-1))
       z          = C._initVars(z,'centers:niveaux_temps',float(niveauTps))
       if niveauTps > time_lvl_max: time_lvl_max = niveauTps
    dt = dt * time_lvl_max
    tc = X.setInterpData2 (t, tc, nature=1, loc='centers', storage='inverse',sameName=1, method='lagrangian',dim=2)
     
t1c= Internal.copyRef(tc)
test.testT(t1c, 5)

t = C.initVars(t, 'centers:cellN', 1.) # init pour les IBCs

# Blanking
s2     = T.addkplane(s)
bodies = [[s2]]
BM     = numpy.array([[1]],numpy.int32)
t      = X.blankCells(t, bodies, BM, blankingType='center_in',dim=2)
t      = X.setHoleInterpolatedPoints(t, depth=+2)
tc     = C.cpVars(t, 'centers:cellN', tc, 'cellN')

sc = C.node2Center(s2)

# Dist2Walls
t = DTW.distance2Walls(t, [sc], type='ortho', loc='centers', signed=1,dim=2)
t = C.center2Node(t, 'centers:TurbulentDistance')

# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
tc = X.setIBCData(t, tc, loc='centers', nature=1, storage='inverse', hi=dh, he=dh,method='lagrangian',dim=2)
t = C.rmVars(t,['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance',\
                'TurbulentDistance','centers:TurbulentDistance'])
tc = C.rmVars(tc, 'cellN') # tres important pour l'instant

t1c= Internal.copyRef(tc)
test.testT(t1c, 6)

##Local time stepping
list_names = ['PointRange','PointRangeDonor','DirDonneur','DirReceveur','Transform','PointPivot','Profondeur','LevelZRcv','LevelZDnr', 'NMratio','DnrZoneName']
if explicit_select == "explicit_local":
    zones = Internal.getNodesFromType(tc, 'Zone_t')
    for z in zones:
        zones2 = Internal.getNodesFromType(z, 'ZoneSubRegion_t')
        for z2 in zones2:
            if z2[0][0:3]=='IBC':
                node = Internal.createNode('PointRange', 'IndexArray_t', value=[0,0,0,0,0,0], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('PointRangeDonor', 'IndexArray_t', value=[0,0,0,0,0,0], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('DirDonneur', 'IndexArray_t', value=[0], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('DirReceveur', 'IndexArray_t', value=[0], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('Transform', 'IndexArray_t', value=[1,2,3], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('PointPivot', 'IndexArray_t', value=[0,0,0], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('Profondeur', 'IndexArray_t', value=[2], children=[])
                Internal.addChild(z2, node, pos=-1)                
                node = Internal.createNode('NMratio', 'IndexArray_t', value=[1,1,1], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('DnrZoneName', 'IndexArray_t', value=[Internal.getValue(z2)], children=[])
                Internal.addChild(z2, node, pos=-1)
                
                t_zone     = Internal.getNodesFromName(t, Internal.getValue(z2))
                time_level = C.getMaxValue(t_zone, 'centers:niveaux_temps')

                node = Internal.createNode('LevelZRcv', 'IndexArray_t', value=[time_level], children=[])
                Internal.addChild(z2, node, pos=-1)
                node = Internal.createNode('LevelZDnr', 'IndexArray_t', value=[time_level], children=[])
                Internal.addChild(z2, node, pos=-1)
# Solver settings
numb = {}
numz = {}
numb["omp_mode"]           = 0
numb["modulo_verif"]       = 10
numz["scheme"]             = "roe_min"
numz["time_step"]          = dt
if explicit_select == 'explicit_local':
    numb["temporal_scheme"]= "explicit_local"
else:
    numb["temporal_scheme"]= "explicit"   

FastC._setNum2Zones(t, numz)
FastC._setNum2Base(t, numb)
    
#C.convertPyTree2File(tc, "tc.cgns")
#C.convertPyTree2File(t ,"mesh.cgns")

(t, tc, metrics) = FastS.warmup(t, tc)

t1= Internal.copyRef(t)
t1c= Internal.copyRef(tc)
Internal._rmNodesFromName(t1, 'Parameter_int')
Internal._rmNodesFromName(t1, 'Parameter_real')
Internal._rmNodesFromName(t1, '.Solver#dtloc')
test.testT(t1, 2)
test.testT(t1c, 3)
#C.convertPyTree2File(t, "Postout.cgns")

nit = 100
for it in range(nit):
    FastS._compute(t, metrics, nit, tc,layer="Python")
    if it%100 == 0:
        print('- %d -'%it); 
        FastS.displayTemporalCriteria(t, metrics, it)

Internal._rmNodesFromName(t, 'Parameter_int')
Internal._rmNodesFromName(t, 'Parameter_real')
Internal._rmNodesFromName(t, '.Solver#dtloc')
test.testT(t, 1)
#C.convertPyTree2File(t, "out.cgns")
