# - TGV FastS (pyTree) -
# - Explicit Local -
#-------------------------
import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Initiator.Adim as Adim
import Generator.PyTree as G
import KCore.test as test
import math

##Please refore to the non-regression test case "explicitLocal_dxdtchange.py" for a more detailed explanation of some settings
##The same approach/logic is used for the current non-regression test case

explicit_select = 'explicit_local'
isPrintAnalysis = False
T0    = 294.
rho0  = 1.2
M     = 0.1
U0    = 34.
Re    = 1600.
L0    = 7e-04
KA    = rho0*U0*U0/16.
R     = 287.053
Mus   = (rho0*U0*L0)/Re
P0    = rho0*R*T0
tfinal = 14.*L0/U0/2.
dt     = 2.*0.0006*L0/U0

if explicit_select == 'explicit_local':
    dt     = 2*dt   
nit    = 100 #int(tfinal/dt)

N      = 66 #Chosen to be a multiple of 3 for the three domains
dx     = 2*math.pi*L0/N
M      = int(N/2*2/3)

if explicit_select == 'explicit_local':
    #### Maillage ####
    a1 = G.cart((-math.pi*L0   ,-math.pi*L0,-math.pi*L0),(dx  ,dx,dx), (  M+1+4  ,N+1,N+1))
    a2 = G.cart((-math.pi*L0/3.,-math.pi*L0,-math.pi*L0),(dx/2,dx,dx), (2*M+1+2*4,N+1,N+1))
    a3 = G.cart(( math.pi*L0/3.,-math.pi*L0,-math.pi*L0),(dx  ,dx,dx), (  M+1    ,N+1,N+1))
    
    a1 = C.addVars(a1, 'centers:niveaux_temps')
    a1 = C.initVars(a1, 'centers:niveaux_temps',1)
    a1 = C.addBC2Zone(a1, 'ovmax', 'BCOverlap', 'imax')
    
    a2 = C.addVars(a2, 'centers:niveaux_temps')
    a2 = C.initVars(a2, 'centers:niveaux_temps',2)
    a2 = C.addBC2Zone(a2, 'ovmin', 'BCOverlap', 'imin')
    a2 = C.addBC2Zone(a2, 'ovmax', 'BCOverlap', 'imax')

    a3 = C.addVars(a3, 'centers:niveaux_temps')
    a3 = C.initVars(a3, 'centers:niveaux_temps',1)
    a3 = C.addBC2Zone(a3, 'ovmin', 'BCOverlap', 'imin')
    ### Arbre t ####
    t = C.newPyTree(['Base']) ; t[2][1][2] += [a1,a2,a3]
    t = X.applyBCOverlaps(t, depth=2)
else:
    #### Maillage ####
    a1 = G.cart((-math.pi*L0   ,-math.pi*L0,-math.pi*L0),(dx  ,dx,dx), (  M+1,N+1,N+1))
    a2 = G.cart((-math.pi*L0/3.,-math.pi*L0,-math.pi*L0),(dx/2,dx,dx), (2*M+1,N+1,N+1))
    a3 = G.cart(( math.pi*L0/3.,-math.pi*L0,-math.pi*L0),(dx  ,dx,dx), (  M+1,N+1,N+1))
    ### Arbre t ####
    t = C.newPyTree(['Base']) ; t[2][1][2] += [a1,a2,a3]    

#-------------------------
# State 
#-------------------------
t = C.addState(t, 'GoverningEquations', 'NSLaminar')
t = C.addState(t, 'EquationDimension' , 3)
t = C.addState(t, adim='dim4',UInf=U0, TInf=T0, PInf=P0, LInf=L0, Mus=Mus)
state      = Adim.dim4(UInf=U0, TInf=T0, PInf=P0, LInf=L0, Mus=Mus)

print("Checking state...start")
list_state = ["RoInf", "RouInf", "RovInf", "RowInf", "RoEInf", "PInf", "TInf", "cvInf", "MInf", "ReInf", "Cs", "Gamma", "RokInf", "RoomegaInf", "RonutildeInf", "Mus", "Cs", "Ts, Pr"]
for i in range(0,len(list_state)):
    print(list_state[i],state[i])
print("Checking state...end")

#-------------------------
# Initialization
#-------------------------
C._initVars(t, '{centers:Density}=(%g+%g*(cos(2.*{centers:CoordinateX}/%g)+cos(2.*{centers:CoordinateY}/%g))*(cos(2.*{centers:CoordinateZ}/%g)+2.))/%g'%(P0,KA,L0,L0,L0,R*T0))
C._initVars(t, '{centers:VelocityX}= %g*sin({centers:CoordinateX}/%g)*cos({centers:CoordinateY}/%g)*cos({centers:CoordinateZ}/%g)'%(U0,L0,L0,L0))
C._initVars(t, '{centers:VelocityY}= -%g*cos({centers:CoordinateX}/%g)*sin({centers:CoordinateY}/%g)*cos({centers:CoordinateZ}/%g)'%(U0,L0,L0,L0))
C._initVars(t, '{centers:VelocityZ} = 0.')
C._initVars(t, '{centers:Temperature} = %g'%T0)
  
if explicit_select != 'explicit_local':
    t = X.connectMatch(t, tol=1.e-6, dim=3)
dist = 2*math.pi*L0 

#-------------------------
# Connections - Raccords
#-------------------------
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],translation=[dist , 0.  ,   0.],dim=3)
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],translation=[   0., dist,   0.],dim=3)
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],translation=[   0., 0.  ,dist ],dim=3)

C.addState2Node__(t, 'EquationDimension', 3)
t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
tc = C.node2Center(t)
if explicit_select == 'explicit_local':
    tc = X.setInterpData2(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, method='lagrangian',dim=3)
else:
    tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, method='lagrangian',dim=3)
tc = C.rmVars(tc, 'FlowSolution')
tc = C.rmVars(tc, 'CellN')

#-------------------------
# Compute
#-------------------------
# Numerics
numb = {}
numz = {}
numb["omp_mode"]           = 0
numb["modulo_verif"]       = 10
numz["scheme"]             = "senseur"
numz["time_step"]          = dt
if explicit_select == 'explicit_local':
    numb["temporal_scheme"]= "explicit_local"
    numb["rk"]        	   = 3
    numb["exp_local"]	   = 2 
    numz["niveaux_temps"]  = 2 # Explicit local : nbre niveaux en temps (=1 si explicit global)
else:
    numb["temporal_scheme"]= "explicit"   
    numz["niveaux_temps"]  = 1 # Explicit local : nbre niveaux en temps (=1 si explicit global)

Fast._setNum2Zones(t, numz);
Fast._setNum2Base(t, numb)
if explicit_select == 'explicit_local':
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    for z in zones:   
        solcenter = Internal.getNodeFromName1(z, 'FlowSolution#Centers') 
        niveau = Internal.getNodeFromName(solcenter, 'niveaux_temps')[1][0][0][0]
        dtloc = Internal.getNodeFromName1(z, '.Solver#define')  # noeud
        level = Internal.getNodeFromName1(dtloc, 'niveaux_temps')  # noeud
        Internal.setValue(level,int(niveau))
            
(t, tc, metrics)  = FastS.warmup(t, tc)

t1 = Internal.copyRef(t)
Internal._rmNodesByName(t1, '.Solver#Param')
Internal._rmNodesByName(t1, '.Solver#ownData')
test.testT(t1, 2)
t1c = Internal.copyRef(tc)
Internal._rmNodesByName(t1c, '.Solver#Param')
Internal._rmNodesByName(t1c, '.Solver#ownData')
test.testT(t1, 3)

times= 0.
timeStep = numz['time_step']

if isPrintAnalysis:
    FileEnstInt=open('enstrophie_N%s_%s.dat'%(N,explicit_select),'a')
    if explicit_select == 'explicit_local':
        C.convertPyTree2File(t,"postwarmup_tgv_explicit_local.cgns")
    else:
        C.convertPyTree2File(t,"postwarmup_tgv_explicit_global.cgns")
        
for it in range(nit):
    FastS._compute(t, metrics, it, tc,layer='Python')
    times += timeStep
    if it%100 == 0:
        print('- %d - %f'%(it, times))
        if isPrintAnalysis:
            (enstrophie, tke) = FastS._computeEnstrophy(t, metrics,times)
            tke = tke/(U0**2)/rho0
            enstrophie = L0**2*enstrophie/U0**2/rho0
            FileEnstInt.write('%s %s %s\n'%(times*U0/L0,tke, enstrophie))
            FileEnstInt.flush()
        if it%1000 == 0:
            FastS.display_temporal_criteria(t, metrics, it)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

test.testT(t, 1)
#C.convertPyTree2File(t, "out.cgns")
#if isPrintAnalysis:
#    if explicit_select == 'explicit_local':
#        C.convertPyTree2File(t,"restart_tgv_explicit_local.cgns")
#    else:
#        C.convertPyTree2File(t,"restart_tgv_explicit_global.cgns")
