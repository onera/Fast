# - compute (pyTree) -
# - Monobloc laminar flat plate -
import Converter.PyTree as C
import FastC.PyTree as Fast
import FastS.PyTree as FastS
import Connector.PyTree as X
import Converter.Internal as Internal
import Transform.PyTree as T
import KCore.test as test

a = C.convertFile2PyTree('Pplane_136_96.cgns')

t = T.subzone(a, (1,20,1), (-1,-1,-1))
t = C.addState(t, MInf=0.2, ReInf=25.e6, MutSMuInf=15)

#recuperation BC disparu par subzone
bc1 = Internal.getNodeFromName(a,'BCWall.10')
Internal.setValue(bc1,'BCWallModel')
Internal.getNodeFromName(t,'ZoneBC')[2].append(bc1)
bc2 = Internal.getNodeFromName(a,'BCFarfieldFill3')
Internal.getNodeFromName(t,'ZoneBC')[2].append(bc2)

coordy = Internal.getNodeFromName(t,'CoordinateY')[1]
dist = Internal.getNodeFromName(t,'TurbulentDistance')[1]
dist[:,:,:]-=coordy[0,2,0]
coordy[:,:,:]-=coordy[0,2,0]
coordy[:,1,:]=coordy[:,2,:]
coordy[:,0,:]=coordy[:,2,:]

for tree in [t]:
   Internal._rmNodesByName(tree, 'ReferenceState')
   Internal._rmNodesByName(tree, 'GoverningEquations')
   for b in Internal.getBases(tree):
     C._addState(b, 'EquationDimension', 2)
     C._addState(b, 'GoverningEquations', 'NSTurbulent')
     C._addState(b, UInf=100, RoInf=1.225, PInf=101325.64, LInf=1, alphaZ=0., adim='dim3')
     for z in Internal.getZones(b):
       sol = Internal.getNodeFromName1(z,'FlowSolution#Centers')
       C._initVars(z,"{centers:ViscosityEddy}=0.000001")
       for sub in ['','_M1','_P1']:
         C._initVars(z,"{centers:Density"+sub+"}=1.225")
         C._initVars(z,"{centers:VelocityX"+sub+"}=100.0")
         C._initVars(z,"{centers:VelocityY"+sub+"}=0")
         C._initVars(z,"{centers:VelocityZ"+sub+"}=0")
         C._initVars(z,"{centers:Temperature"+sub+"}=288.05")
         C._initVars(z,"{centers:TurbulentSANuTilde"+sub+"}=0.000001")

#C.convertPyTree2File(t, 'verif.cgns')
#stop

sample=3
modulo_verif=100
numb = {'temporal_scheme':'implicit', 'ss_iteration':1, 'modulo_verif':modulo_verif }
numz = {'ssdom_IJK': [600, 30, 20000], 
       'io_thread': -4, 
       'time_step': 0.00000152280, 
       'time_step_nature': 'local', 
       'cfl': 100, 
       'scheme':'roe', 
       'slope':'o3', 
       'wallmodel_sample':sample,
       'epsi_newton': 0.02}
Fast._setNum2Zones(t, numz) ; Fast._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

t1=Internal.copyTree(t)
Internal._rmNodesByName(t1, '.Solver#Param')
Internal._rmNodesByName(t1, '.Solver#ownData')
test.testT(t1, 1)
#stop

nit = 500; time = 0.
#nit = 1; time = 0.
for it in range(nit):
    FastS._compute(t, metrics, it)
    if it%modulo_verif == 0:
        print('- %d - %g'%(it, time))
        FastS.display_temporal_criteria(t, metrics, it)
    time += numz['time_step']

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 2)
