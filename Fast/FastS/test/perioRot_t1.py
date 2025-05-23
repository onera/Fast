# - cylinder (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Converter.Internal as Internal
import KCore.test as test

angle = 45.
a = G.cylinder((0.,0.,0.), 0.5, 1., angle , 0., 2., (50,50,5))
t = C.newPyTree(['Base',a])
C._addState(t, 'GoverningEquations', 'NSLaminar')
C._addState(t, MInf=0.1, ReInf=1600., adim='adim2funk')

# Get dim
dim = 3
node = Internal.getNodeFromName(t, 'EquationDimension')
if node is not None: dim = Internal.getValue(node)

# Solution initiale
eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
Model = 'NSLaminar'
if eqs is not None: Model = Internal.getValue(eqs)

ret = C.isNamePresent(t, 'centers:Density')
if ret != 1: # Density not present
    state = Internal.getNodeFromType(t, 'ReferenceState_t')
    if state is None:
        raise ValueError('Reference state is missing in input cgns.')
    vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
            'EnergyStagnationDensity']
    for v in vars:
        node = Internal.getNodeFromName(state, v)
        if node is not None:
            val = float(node[1][0])
            C._initVars(t, 'centers:'+v, val)
        else:
            raise ValueError(v + ' is missing in ReferenceState.')
    if Model == 'NSTurbulent':
        vars = ['TurbulentSANuTildeDensity']
        for v in vars:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                C._initVars(t, 'centers:'+v, val)
ret = C.isNamePresent(t, 'centers:TurbulentDistance')
if Model == 'NSTurbulent' and ret != 1: # Wall distance not present
    import Dist2Walls.PyTree as DTW
    walls = C.extractBCOfType(t, 'BCWall')
    DTW._distance2Walls(t, walls, loc='centers', type='ortho')

t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,angle])
C._fillEmptyBCWith(t, 'Extrap', 'BCExtrapolate', dim=3)
#t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period kmin', 'BCautoperiod', 'kmin')
#t[2][1][2][0] = C.addBC2Zone(t[2][1][2][0], 'period kmax', 'BCautoperiod', 'kmax')


# Ajout des ghost cells
C.addState2Node__(t, 'EquationDimension', dim)
NGhostCells = 2
Internal._addGhostCells(t, t, NGhostCells, adaptBCs=1, fillCorner=0)
if dim == 2:
    T._addkplane(t)
    T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.01)
    T._makeDirect(t)
#C.convertPyTree2File(t,"in.cgns")
#C._initVars(t,"centers:cellN",1.)

tc = C.node2Center(t)
tc = X.setInterpData(t,tc, storage='inverse',loc='centers',penalty=1,nature=1,itype='abutting')
C._rmVars(tc, 'FlowSolution')

# initialisation
C._initVars(t, '{centers:Rayon}=sqrt({centers:CoordinateX}**2+{centers:CoordinateY}**2) ')
C._initVars(t, '{centers:Density}= 1.')
C._initVars(t, '{centers:VelocityX}= {centers:CoordinateX}/{centers:Rayon}')
C._initVars(t, '{centers:VelocityY}= {centers:CoordinateY}/{centers:Rayon}')
C._initVars(t, '{centers:VelocityZ} = 1.')
C._initVars(t, '{centers:Temperature}=1.')

# Numerics
numb = {}; numz = {}
numb["temporal_scheme"]    = "explicit"
FastC._setNum2Base(t, numb); FastC._setNum2Zones(t, numz)

(t, tc, metrics) = FastS.warmup(t, tc, graph=None)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t,1)

#C.convertPyTree2File(t, 'out.cgns')
