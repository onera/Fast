# - cylinder (pyTree) -
import Transform.PyTree as T
import Initiator.PyTree as I
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
import math


#grille cylindrique d'axe Z
a = G.cart((0,0,0), (0.5,0.5,0.25), (200,100,15))

t = C.newPyTree(['Base',a])

U0= 140.67
P0= 89280.81
L0=1.  
R0=1.111711
T0=279.15
t = C.addState(t, 'GoverningEquations', 'NSLaminar')
t = C.addState(t, UInf=U0, RoInf=R0, PInf=P0, LInf=L0, alphaZ=0., adim='dim3')

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
        raise ValueError, 'Reference state is missing in input cgns.'
    vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
            'EnergyStagnationDensity']
    for v in vars:
        node = Internal.getNodeFromName(state, v)
        if node is not None:
            val = float(node[1][0])
            C._initVars(t, 'centers:'+v, val)
        else:
            raise ValueError, v + ' is missing in ReferenceState.'



# Ajout des ghost cells
C.addState2Node__(t, 'EquationDimension', dim)
NGhostCells = 2
t = Internal.addGhostCells(t, t, NGhostCells, adaptBCs=1, fillCorner=0)

#import sys
#sys.exit()



# initialisation

t= C.initVars(t, '{centers:Density}= 1.')
t= C.initVars(t, '{centers:VelocityX}= 1.')
t= C.initVars(t, '{centers:VelocityY}= 1.')
t= C.initVars(t, '{centers:VelocityZ}= 1.')
t= C.initVars(t, '{centers:Temperature}=1.')

zones = Internal.getZones(t)

for z in  zones:

   xx = Internal.getNodeFromName2( z, "CoordinateX" )[1]
   yy = Internal.getNodeFromName2( z, "CoordinateY" )[1]
   zz = Internal.getNodeFromName2( z, "CoordinateZ" )[1]

   sol = Internal.getNodeFromName1( z, "FlowSolution#Centers" )
   vx = Internal.getNodeFromName1( sol, "VelocityX" )[1]
   vy = Internal.getNodeFromName1( sol, "VelocityY" )[1]
   vz = Internal.getNodeFromName1( sol, "VelocityZ" )[1]
   
   #for k in range( xx.shape[2] ):
   #  xx[10,10,k]=  xx[10,10,k]+0.000001

   for k in range( vx.shape[2] ):
     for j in range( vx.shape[1] ):
       for i in range( vx.shape[0] ):
          x1 =  xx[i,j,k]
          y1 =  yy[i,j,k]
          z1 =  zz[i,j,k]

          vx[i,j,k]= 10.*math.cos(x1)*math.sin(y1)*math.sin(z1)
          vy[i,j,k]= 10.*math.sin(x1)*math.cos(y1)*math.sin(z1)
          vz[i,j,k]= 10.*math.cos(x1)*math.sin(y1)*math.cos(z1)


'''   k=21
   for j in range( 4 ):
       for i in range( 4 ):
          x1 =  xx[i,j,k]
          y1 =  yy[i,j,k]
          z1 =  zz[i,j,k]
          print 'VV', vx[i,j,k],vy[i,j,k],vz[i,j,k],i,j,k
          print 'XX', math.cos(x1),math.sin(y1),math.sin(z1)

import sys; sys.exit()
'''
# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numz = {}
numz["time_step_nature"]   = 'global'
numz["time_step"]          = 0.8e-9
numz["sgsmodel"]           = 'smsm'
numz["io_thread"]          =  0
numz["scheme"]             = "senseur" #ausmpred ou senseur
#numz["cache_blocking_I"]   = 2000
#numz["cache_blocking_J"]   = 400
#numz["cache_blocking_K"]   = 400

Fast._setNum2Zones(t, numz)
Fast._setNum2Base(t, numb)

tc=None
(t, tc, metrics) = FastS.warmup(t, tc, graph=None)

#t = Internal.rmGhostCells(t, t, 1 )
test.testT(t, 1)
import sys;sys.exit()

ref = C.convertFile2PyTree( "Data/sgs_t1.ref1")

zones_ref = Internal.getZones(ref)
ro_ref = Internal.getNodeFromName2( zones_ref[0], "ViscosityEddy" )[1]
zones = Internal.getZones(t)
ro = Internal.getNodeFromName2( zones[0], "ViscosityEddy" )[1]

print 'shape ro', ro.shape

'''   for j in range( ro.shape[1] ):
    for i in range( ro.shape[0] ):
      if( abs(ro_ref[i,j,k] - ro[i,j,k])) > 1.e-7:
'''
i=10
j=82
#for j in range( 78, ro.shape[1]-10 ):
for k in range( ro.shape[2] ):
             #print i-1,j-1,k-1, ro_ref[i,j,k] - ro[i,j,k]
             print k-1, j-1, ro[i,j,k]
      
             #ro[i,j,k] = ro[i,j,k] - ro_ref[i,j,k]
C.convertPyTree2File(t, "bug.cgns")
import sys;sys.exit()



nit = 1; time = 0.
timeStep = numz['time_step']
for it in xrange(nit):
    FastS._compute(t, metrics, it)
    time += timeStep


C.convertPyTree2File(t, 'out.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

