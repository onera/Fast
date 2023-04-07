# - cylindre chimere avec mouvement -
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import RigidMotion.PyTree as R
import Geom.PyTree as D
import Post.PyTree as P
import Connector.ToolboxIBM as TIBM
import Generator.PyTree as G
import Distributor2.PyTree as D2
import Transform.PyTree as T
import Compressor.PyTree as Compressor
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Connector.Mpi as Xmpi
import Converter.Mpi as Cmpi
import CPlot.PyTree as CPlot
import KCore.test as test
import numpy

# case
N = 100; h = 1./(N-1)
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (2*N,int(0.5*N),1)) 
a = X.connectMatch(a, dim=2)
C._addBC2Zone(a, 'wall', 'BCWall', 'jmin')
C._fillEmptyBCWith(a, 'overlap', 'BCOverlap')

N = 280; h = 8./(N-1)
b = G.cart((-5,-4,0), (h,h+0.01,1), (N,N,1))
C._fillEmptyBCWith(b, 'farfield', 'BCFarfield')

t = C.newPyTree(['CYL',a,'FOND',b])

# prepare
for name in ['CYL','FOND']:
    C._addState(Internal.getNodeFromName1(t, name),
                adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=1.e6, 
                EquationDimension=2, GoverningEquations='NSLaminar')

R._setPrescribedMotion3(Internal.getNodeFromName1(t, 'CYL'),
                'rot', axis_pnt=(-1.5,0.,0.), axis_vct=(0,0,1), omega=0.5)

bases = Internal.getBases(t)
tb = C.newPyTree()
for c, b in enumerate(bases):
    bd = C.extractBCOfType(b, 'BCOverlap', shift=15)
    bd = C.convertArray2Tetra(bd)
    if bd != []:
        bd = T.join(bd)
        bd[0] = 'BODY%d'%c
        n = Internal.getNodeFromName(b, 'TimeMotion')
        if n is not None: bd[2].append(n)
        bb = Internal.copyRef(b)
        bb[2] = [bd]
        tb[2].append(bb)
R._copyGrid2GridInit(tb)

Internal._addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
t = T.addkplane(t)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
T._makeDirect(t)

tc = C.node2Center(t)
zones = Internal.getZones(Internal.getNodeFromName1(t, 'CYL'))
zonesC = Internal.getZones(Internal.getNodeFromName1(tc, 'CYL'))
X._setInterpData(zones, zonesC, nature=1, loc='centers', storage='inverse', sameName=1, dim=2)
zones = Internal.getZones(Internal.getNodeFromName1(t, 'FOND'))
zonesC = Internal.getZones(Internal.getNodeFromName1(tc, 'FOND'))
X._setInterpData(zones, zonesC, nature=1, loc='centers', storage='inverse', sameName=1, dim=2)
C._rmVars(tc, 'FlowSolution')

C._initVars(t, 'centers:cellN=1.')
vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity']

for b in Internal.getBases(t):
    state = Internal.getNodeFromType(b, 'ReferenceState_t')
    Model = Internal.getNodeFromName(b, 'GoverningEquations')
    Model = Internal.getValue(Model)
    if Model == 'NSTurbulent': allvars = vars + ['TurbulentSANuTildeDensity']
    else: allvars = vars
    for v in allvars:
        node = Internal.getNodeFromName(state, v)
        if node is not None:
            val = float(node[1][0])
            C._initVars(b, 'centers:'+v, val)
            
for b in Internal.getBases(t):
    Model = Internal.getNodeFromName(b, 'GoverningEquations')
    Model = Internal.getValue(Model)
    if Model == 'NSTurbulent': # Wall distance
        import Dist2Walls.PyTree as DTW
        walls = C.extractBCOfType(b, 'BCWall')
        if walls != []:
            DTW._distance2Walls(b, walls, loc='centers', type='ortho')
        else: C._initVars(b, 'centers:TurbulentDistance', 1000000.)

R._copyGrid2GridInit(t)
R._copyGrid2GridInit(tc)

# compute
numb={"temporal_scheme": "explicit", "omp_mode":0, "modulo_verif":10}
numz={"time_step": 2.e-3, "scheme":"roe_min"}

Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)

zones = Internal.getZones(t)
for z in zones:
    timeMotion = Internal.getNodeFromName(z, 'TimeMotion')
    if timeMotion is not None:
        define = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(define, 'motion', 'DataArray_t', value='deformation')

(t, tc, metrics) = FastS.warmup(t, tc)

C._initVars(t, "{centers:cellN#Motion}=1.")
X._applyBCOverlaps(Internal.getNodeFromName1(t,'CYL'), depth=4, val=2, cellNName='cellN#Motion')
X._applyBCOverlaps(Internal.getNodeFromName1(t,'CYL'), depth=2, val=0, cellNName='cellN#Motion')
C._initVars(t, "{centers:cellN#MotionInit}={centers:cellN#Motion}")

nbases = len(Internal.getBases(t))
nbodies = len(Internal.getBases(tb))
BM = numpy.zeros((nbases, nbodies),dtype=numpy.int32)
BM[1,0] = 1.

tBB = Cmpi.createBBoxTree(t)

intersectionDict={}
for z1 in Internal.getZones(Internal.getNodeFromName1(tBB, 'CYL')):
    for z2 in Internal.getZones(Internal.getNodeFromName1(tBB, 'FOND')):
        Fast._addPair(intersectionDict, z1[0], z2[0])

procDict=None; graphX={}
procDict = Cmpi.getProcDict(tBB)
graphX = Cmpi.computeGraph(tBB, type='bbox2', t2=None, 
                           procDict=procDict, reduction=False,
                           intersectionsDict=intersectionDict)

dictOfADT={}
(dictOfNobOfRcvZones,dictOfNozOfRcvZones)   = Fast.getDictOfNobNozOfRcvZones(t, intersectionDict)
(dictOfNobOfRcvZonesC,dictOfNozOfRcvZonesC) = Fast.getDictOfNobNozOfRcvZones(tc, intersectionDict)
(dictOfNobOfDnrZones,dictOfNozOfDnrZones)   = Fast.getDictOfNobNozOfDnrZones(tc, intersectionDict, dictOfADT)

time = 0

for it in range(500):
    time += time_step

    R._evalPosition(tb, time)
    R._evalPosition(t, time)
    R._evalPosition(tc, time)
    R._evalGridSpeed(t, time)
    FastS.copy_velocity_ale(t, metrics, it=it)

    C._initVars(t,"{centers:cellN#Motion}={centers:cellN#MotionInit}")
    bodies=[]
    for base in Internal.getBases(tb):
        bodies.append(Internal.getZones(base))
    X._blankCells(t, bodies, BM, cellNName='cellN#Motion', XRaydim1=100, XRaydim2=100, dim=2)
    X._setHoleInterpolatedPoints(Internal.getNodeFromName1(t, 'FOND'), depth=2, loc='centers', addGC=False, cellNName='cellN#Motion', dir=0)
    C._initVars(t, "{centers:cellN}=minimum(2,{centers:cellN#Motion})")

    ucData = (graphX, intersectionDict, dictOfADT, 
            dictOfNobOfRcvZones, dictOfNozOfRcvZones,
            dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
            dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC,
            time, procDict, True, 0, 1, 2)
    
    FastS._compute(t, metrics, it, tc, None, layer="Python", ucData=ucData)

for adt in dictOfADT.values():
    if adt is not None: C.freeHook(adt)
C._rmVars(t, 'centers:cellN#Motion')
C._rmVars(t, 'centers:cellN#MotionInit')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData') 

test.testT(t, 1)

