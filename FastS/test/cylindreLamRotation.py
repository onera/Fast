# - cylinder (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Converter.Internal as Internal
import Fast.Utils as Utils
import KCore.Adim as Adim
import KCore.test as test
import numpy
import math 
import sys  
import script_prepareCaseHdf_r37 as SPCH

angle = 10.

anglex =angle
angley = 0.
anglez = 0.

a = G.cylinder((0.,0.,0.), 0.24, 0.3, angle , 0., 0.22, (21,77,35)) 

b = T.rotate(a, (0.,0.,0.), (0.,1.,0.), 90.); b[0] = 'cart'

b = T.reorder(b, (2,-3,-1))
t = C.newPyTree(['Base',b])

dist = G.cart((0,0,0), (1./76, 1,1), (77,1,1) )
dist = G.enforceMoinsX(dist, 2.e-4, (10,10))
dist = G.enforcePlusX(dist, 2.e-4, (10,10))

b = G.map( b, dist, dir=3)

t = C.newPyTree(['Base',b])

size_fenbc=1920

dim=3
t = C.addState(t, 'EquationDimension', dim)
t = C.addState(t, 'GoverningEquations', 'NSTurbulent')
#C.addState2Node__(t[2], 'EquationDimension', dim)
#C.addState2Node__(t[2], 'GoverningEquations', 'NSTurbulent')
t = C.addState(t, MInf=0.5, ReInf=4000000., UInf= 166.057, TInf= 274.476, PInf= 85418 , LInf=4.,  adim='dim1')

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


t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.], rotationAngle=[anglex,-angley,anglez])
#BC kmin
C._addBC2Zone(t[2][1][2][0], 'wall_fixe', 'BCWall', 'kmin')
Prop = Internal.createChild(  Internal.getNodesFromName( t[2][1][2][0] ,'wall_fixe')[0]  ,'.Solver#Property','UserDefinedData_t')
Internal.createChild(Prop,'mobile_coef','DataArray_t',value=0.)

#BC kmax
C._addBC2Zone(t[2][1][2][0], 'wall_mobile',   'BCWall', 'kmax')
Prop = Internal.createChild(  Internal.getNodesFromName( t[2][1][2][0] ,'wall_mobile')[0]  ,'.Solver#Property','UserDefinedData_t')
Internal.createChild(Prop,'mobile_coef','DataArray_t',value=1.)

#BC imax
C._addBC2Zone(t[2][1][2][0], 'Inj1',   'BCInj1', 'imax')
Prop = Internal.createChild(  Internal.getNodesFromName( t[2][1][2][0] ,'Inj1')[0]  ,'.Solver#Property','UserDefinedData_t')
inj1 = numpy.empty( size_fenbc, numpy.float64); inj1[0:size_fenbc] = 1.
inj2 = numpy.empty( size_fenbc, numpy.float64); inj2[0:size_fenbc] = 0.
inj3 = numpy.empty( size_fenbc, numpy.float64); inj3[0:size_fenbc] = 0.
inj4 = numpy.empty( size_fenbc, numpy.float64); inj4[0:size_fenbc] = 101325.
inj5 = numpy.empty( size_fenbc, numpy.float64); inj5[0:size_fenbc] = 289537.
inj6 = numpy.empty( size_fenbc, numpy.float64); inj6[0:size_fenbc] = 3.4e-6
Internal.createChild(Prop,'txv'                ,'DataArray_t', inj1)
Internal.createChild(Prop,'tyv'                ,'DataArray_t', inj2)
Internal.createChild(Prop,'tzv'                ,'DataArray_t', inj3)
Internal.createChild(Prop,'stagnation_pressure','DataArray_t', inj4)
Internal.createChild(Prop,'stagnation_enthalpy','DataArray_t', inj5)
Internal.createChild(Prop,'inj_tur1'           ,'DataArray_t', inj6)

#BC imin
C._addBC2Zone(t[2][1][2][0], 'outpres',   'BCOutpres', 'imin')
Prop = Internal.createChild(  Internal.getNodesFromName( t[2][1][2][0] ,'outpres')[0]  ,'.Solver#Property','UserDefinedData_t')
out1 = numpy.empty( size_fenbc, numpy.float64); out1[0:size_fenbc] = 91192.5
Internal.createChild(Prop,'pressure','DataArray_t', out1)
Internal.createChild(Prop,'k_piv'   ,'DataArray_t',value=1.)

ret = C.isNamePresent(t, 'centers:TurbulentDistance')
if Model == 'NSTurbulent' and ret != 1: # Wall distance not present
    import Dist2Walls.PyTree as DTW
    walls = C.extractBCOfType(t, 'BCWall')
    t = DTW.distance2Walls(t, walls, loc='centers', type='ortho')

# Ajout des ghost cells
NGhostCells = 2
t = Internal.addGhostCells(t, t, NGhostCells, adaptBCs=1, fillCorner=0)

tc = C.node2Center(t)
C.addState2Node__(tc, 'EquationDimension', dim)
C.addState2Node__(tc, 'GoverningEquations', 'NSTurbulent')
tc = X.setInterpData(t,tc, storage='inverse',loc='centers',penalty=1,nature=1,itype='abutting')
C._rmVars(tc, 'FlowSolution')

# initialisation
#C._initVars(t, '{centers:Rayon}=sqrt ({centers:CoordinateY}**2+{centers:CoordinateZ}**2) ')
C._initVars(t, '{centers:Density}= 1.08419')
C._initVars(t, '{centers:VelocityX} = 166.057')
C._initVars(t, '{centers:VelocityY}= 0.')
C._initVars(t, '{centers:VelocityZ}= 0.')
C._initVars(t, '{centers:Temperature}=274.476')




#=================
restart = 'no'
steady  = 'yes'
NIT     = 100
#=================

# Rotation axe Ox avec vitess angulaire en rpm
omgrpm = -17188.*math.pi/30
Naube1 = 36

#NIT  = 1
NTOUR    = 1
NIT_TOUR = 360
RED_DT   = 0.10
time_step = RED_DT*2.*math.pi/(NIT_TOUR*abs(omgrpm))
print'+++++++++++++++++++++++++'
print'time_step = ',time_step
print'niter     e= ',NIT
print'+++++++++++++++++++++++++'

teta      = 0.
tetap     = 0.
time_ale  = 0.

# Auto-restart
it0 = 1
time0 = 0.
it = 0

modulo_verif = 100

# Numerics global
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["modulo_verif"]       =  modulo_verif
numb["ss_iteration"]       =  3
# Numerics par zone
numz = {}
numz["time_step"]          = time_step
numz["motion"]             = "rigid" 
numz["rotation"]           =  [ 1.,0.,0., 0., 0., 0., 0., 0.]  # axe rotation, centre rotation, freq+ampli  
numz["epsi_newton"]        = 0.20
numz["scheme"]             = "ausmpred" 
numz["time_step_nature"]   = "local"
numz["cfl"]                = 10

Fast._setNum2Zones(t, numz)
Fast._setNum2Base(t, numb)


# Number or records to store residuals 
nrec = NIT/modulo_verif

time = time0

if steady == 'no' :
    time_step = Internal.getNodeFromName(t, 'time_step')
    time_step = Internal.getValue(time_step)

if restart == 'yes' : 
    NIT_REP = Internal.getNodeFromName(t, 'Iteration')
    NIT_REP = Internal.getValue(NIT_REP)
else :
    NIT_REP  = 1

nit = NIT+NIT_REP

if steady == 'no' :
	first = Internal.getNodeFromName1(t, 'Time')
	if (first != None): time0 = Internal.getValue(first)
	else:               time0 = 0

time = time0


print'it_init, it_fin :',NIT_REP-1,NIT_REP+NIT-1

#Initialisation parametre vitesse d entrainememnt
time_ale = time + 0.5*time_step
teta = omgrpm*time_ale
tetap= omgrpm
infos_ale=[ teta, tetap]

(t, tc, metrics) = FastS.warmup(t, tc, None, infos_ale)

#C.convertPyTree2File(tc, 'tc_test.cgns')

# Preparation extraction debits 
# => creation arbres t_debit_in et t_debit_out
t_debit_in    = FastS.createStressNodes(t, BC=['BCInj1'])
t_debit_out   = FastS.createStressNodes(t, BC=['BCOutpres'])
t_debit_wall  = FastS.createStressNodes(t, BC=['BCWall'])
#t_debit_perio = FastS.createStressNodes(t, BC=['BCPeriodic'])

node_debit_in   = SPCH.prepareDebit2Fast(t_debit_in)
node_debit_out  = SPCH.prepareDebit2Fast(t_debit_out)
node_debit_wall = SPCH.prepareDebit2Fast(t_debit_wall)
#node_debit_perio= SPCH.prepareDebit2Fast(t_debit_perio)

# Listes iteration/debit amont/debit val
#debit_amont = CV.array('iteration,convflux_ro', nrec, 1, 1)
#debit_aval  = CV.array('iteration,convflux_ro', nrec, 1, 1)
#debit_wall  = CV.array('iteration,convflux_ro', nrec, 1, 1)
#debit_perio = CV.array('iteration,convflux_ro', nrec, 1, 1)


#C.convertPyTree2File(tc, 'tc_test.cgns')

#C.convertPyTree2File(t, 'out.cgns')
#sys.exit()

modulo_debit = 100
modulo_convNewt =  modulo_verif
for it in xrange(NIT_REP, NIT_REP+NIT):
#for it in xrange(NIT):
    if steady == 'no':
	time_ale = time + 0.5*time_step
	teta = omgrpm*time_ale
	tetap= omgrpm
    
	FastS._motionlaw(t, teta, tetap)
    
    FastS._compute(t, metrics, it, tc)
    
    if (it%modulo_convNewt == 0):
      FastS.display_temporal_criteria(t, metrics, it)
    
    

    if (it%modulo_debit == 0):
      # Extraction debits	

      #FastS._applyBC(t,metrics)
      #debit_perio= FastS._computeStress(t, t_debit_perio, metrics)
      #debit_perio= Naube1*SPCH.extractDebit2Fast(node_debit_perio)
      #print  "debit perio =", debit_perio

      debit_in = FastS._computeStress(t, t_debit_in, metrics)
      debit_in = Naube1*SPCH.extractDebit2Fast(node_debit_in)                #we extract the debit by doing a sum of all the cells' contributions 
    
      debit_out = FastS._computeStress(t, t_debit_out, metrics)
      debit_out = Naube1*SPCH.extractDebit2Fast(node_debit_out)
     
      debit_wall = FastS._computeStress(t, t_debit_wall, metrics)
      debit_wall = Naube1*SPCH.extractDebit2Fast(node_debit_wall)
    
    
      #CV.setValue(debit_amont, ((it-NIT_REP+1)/modulo_verif ,1,1), [it,-1*debit_in])
      print  "debit amont =", -1*debit_in
    
      #CV.setValue(debit_aval, ((it-NIT_REP+1)/modulo_verif,1,1), [it,debit_out])
      print  "debit aval =", debit_out
      
      print  "debit wall =", debit_wall, debit_in+debit_out

    
    
    if (it%modulo_debit) == 0:
        print "it,time,omega*t = ",it,time,omgrpm*time*180/math.pi
    
    FastS._movegrid(t)
    
    #~ if (it == NIT-1):
	 #~ FastS._movegrid(t)
	 #~ print "movegrid for post processing"
    
    if steady == 'no' : time += time_step

#sys.exit()


Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=nit)
if steady == 'no' : 
    Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time)

nitm1 = nit-1

#=============================
# Supprime les champs inutiles
#=============================
Internal._rmNodesByName(t, '.Solver#Param')
#Internal._rmNodesByName(t, '.Solver#define')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(t, 'FlowSolution')
#Internal._rmNodesByName(t, 'ZoneBC')
Internal._rmNodesByName(t, 'Iteration')
Internal._rmNodesByName(t, 'GridCoordinates#Init')
#Internal._rmNodesByName(t, 'FlowEquationSet')
#Internal._rmNodesByName(t, 'ReferenceState')
vars = ['centers:Density_M1', 'centers:VelocityX_M1', 'centers:VelocityY_M1', 'centers:VelocityZ_M1', 'centers:Temperature_M1', 'centers:TurbulentSANuTilde_M1''centers:Density_P1', 'centers:VelocityX_P1', 'centers:VelocityY_P1', 'centers:VelocityZ_P1', 'centers:Temperature_P1','centers:TurbulentSANuTilde_P1']
C._rmVars(t, vars)

#C.convertPyTree2File(t, 'out.cgns')

test.testT(t, 1)

