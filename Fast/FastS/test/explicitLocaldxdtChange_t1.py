# - compute (pyTree) -
# - Lamb vortex sur 3 domaines - Explicit Local
# - based on explicitLocal_t1.py but uses domain with different grid size
import Generator.PyTree as G
import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I
import FastC.PyTree as internal
import FastS.PyTree as FastS
import Connector.PyTree as X
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.test as test

##Note: (1)Current non-regression test case is for explicit local
##         However, explicit global (coined explicit henceforth) is also supported to
##         compare the two methods for a single test case
##
##      (2) The domain is greatly increased compard to 'explicitLocal_t1.py' to avoid
##          "ugly" domain boundary effects
##expicit_select = ['explicit_local','explicit']

isPrintAnalysis = False
explicit_select = 'explicit_local'
if explicit_select == 'explicit_local':
    #### Maillage ####
    ## Note: 3 cells is the lower limit but chosen for the curren test case
    ## __|___|___|_i_|i+1|i+2|
    ##               |a|b|c|d|
    ## For interpolation at "a" values at "i" and "i+1" are needed.
    ## therefore a physical overlap of 3 cells is needed for a computational overlap
    ## of 2 cells
    a1 = G.cart((-50,-50,0),(0.4,0.5,0.1), (101+4,201,1))
    a1 = C.addVars(a1, 'centers:niveaux_temps')
    a1 = C.initVars(a1, 'centers:niveaux_temps',1)
    a1 = C.addBC2Zone(a1, 'ovmax', 'BCOverlap', 'imax')

    #Overlap for fine is always twice of coarse
    a2 = G.cart((-10,-50,0),(0.2,0.5,0.1), (101+8,201,1))
    a2 = C.addVars(a2, 'centers:niveaux_temps')
    a2 = C.initVars(a2, 'centers:niveaux_temps',2)
    a2 = C.addBC2Zone(a2, 'ovmin', 'BCOverlap', 'imin')
    a2 = C.addBC2Zone(a2, 'ovmax', 'BCOverlap', 'imax')

    a3 = G.cart((10,-50,0),(0.4,0.5,0.1), (101,201,1))
    a3 = C.addVars(a3, 'centers:niveaux_temps')
    a3 = C.initVars(a3, 'centers:niveaux_temps',1)
    a3 = C.addBC2Zone(a3, 'ovmin', 'BCOverlap', 'imin')

    ### Arbre t ####
    t = C.newPyTree(['Base']) ; t[2][1][2] += [a1,a2,a3]
    ## This sets the computational overlap
    t = X.applyBCOverlaps(t, depth=2)
else:
    #### Maillage ####
    a1 = G.cart((-50,-50,0),(0.4,0.5,0.1), (101,201,1))
    a2 = G.cart((-10,-50,0),(0.2,0.5,0.1), (101,201,1))
    a3 = G.cart(( 10,-50,0),(0.4,0.5,0.1), (101,201,1))
    ### Arbre t ####
    t = C.newPyTree(['Base']) ; t[2][1][2] += [a1,a2,a3]

### Conditions initiales ####
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, 'EquationDimension', 2)
t = C.addState(t, MInf=0.2, ReInf=100.)
t = I.initLamb(t, position=(0.,0.), Gamma=2., MInf=0.2, loc='centers')

### Raccords ####
if explicit_select != 'explicit_local':
    t = X.connectMatch(t, tol=1.e-6, dim=2)
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],translation=[100., 0.,0.],dim=2)
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],translation=[  0.,100,0.],dim=2)


### Ajout des ghost-cells ###
C.addState2Node__(t, 'EquationDimension', 2)
t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
t = T.addkplane(t)
t = T.contract(t, (0,0,0), (1,0,0), (0,1,0),0.025)
t = T.makeDirect(t)

### Construction du tc ###
tc = C.node2Center(t)
if explicit_select == 'explicit_local':
    tc = X.setInterpData3(t, tc, nature=1, loc='centers', storage='inverse',sameName=1, method='lagrangian',dim=2)
else:
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse',sameName=1, method='lagrangian',dim=2)
tc = C.rmVars(tc, 'FlowSolution')
tc = C.rmVars(tc, 'CellN')


# Numerics
dt = 0.02
if explicit_select == 'explicit_local':
    dt = 2*dt
tend= 80
NIT = 1000 #int(tend/dt)

numb = {}
numz = {}
numb["omp_mode"]           = 0
numb["modulo_verif"]       = 10
numz["scheme"]             = "ausmpred"
numz["time_step"]          = dt
if explicit_select == 'explicit_local':
    numb["temporal_scheme"]= "explicit_local"
    numb["rk"]        	   = 3
    numb["exp_local"]	   = 2 
    numz["niveaux_temps"]  = 2 # Explicit local : nbre niveaux en temps (=1 si explicit global)
else:
    numb["temporal_scheme"]= "explicit"   
    numz["niveaux_temps"]  = 1 # Explicit local : nbre niveaux en temps (=1 si explicit global)

internal._setNum2Zones(t, numz) 
internal._setNum2Base(t, numb)

if explicit_select == 'explicit_local':
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    for z in zones:   
        solcenter = Internal.getNodeFromName1(z, 'FlowSolution#Centers') 
        niveau = Internal.getNodeFromName(solcenter, 'niveaux_temps')[1][0][0][0]
        dtloc = Internal.getNodeFromName1(z, '.Solver#define')  # noeud
        level = Internal.getNodeFromName1(dtloc, 'niveaux_temps')  # noeud
        Internal.setValue(level,int(niveau))

(t,tc,metrics) = FastS.warmup(t,tc)

nit=NIT; times = 0.
timeStep = numz['time_step']
for it in range(nit):
    FastS._compute(t, metrics, it,tc,layer='Python')
    if it%100 == 0:
        print('- %d - %f'%(it, times))
        FastS.display_temporal_criteria(t, metrics, it)
    times += timeStep

#t = P.computeGrad2(t, 'centers:Density')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

test.testT(t, 1)
#C.convertPyTree2File(t,"out.cgns")

##Note: Uncomment below to test explicit global vs explicit local
#if explicit_select == 'explicit_local':
#    C.convertPyTree2File(t,"restart_explicit_local.cgns")
##    #Uncomment below to compare curve in Visit
##    t=FastS.rmGhostCells(t,2)
##    C.convertPyTree2File(t,"restart_explicit_local.plt")
#else:
#    C.convertPyTree2File(t,"restart_explicit_global.cgns")
##    #Uncomment below to compare curve in Visit
##    t=FastS.rmGhostCells(t,2)
##    C.convertPyTree2File(t,"restart_explicit_global.plt")

