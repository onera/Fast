''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''                   CONVECTED VORTEX NS/LBM 2 DOMAINS                      '''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import FastASLBM.PyTree as FastLBM
import Connector.PyTree as X
import Post.PyTree as P
import Initiator.PyTree as I
import Transform.PyTree as T
import Apps.Fast.Couplage_LBMNS as App

import KCore.test as test

import math
import numpy
import time
import matplotlib.pyplot as plt

myApp = App.CLBMNS(format='single')

#==========================================
# Generation du maillage
#==========================================
x_min = 0.0
L = 1.
Nx=101; Ny = Nx; Nz = 5
n = Nx-1
dx = L/n
z_min = -4*dx
a1 = G.cart((x_min,x_min,z_min), (dx,dx,dx), (Nx,Ny,Nz))
a2 = G.cart((x_min+L,x_min,z_min), (dx,dx,dx), (Nx,Ny,Nz))
t = C.newPyTree(['Base',a1,a2])

SOLVER_TYPE = [1,0]
t,tc = myApp.prepare(t, t_out=None, tc_out=None, solvertype=SOLVER_TYPE, translation=[2*(Nx-1)*dx,(Ny-1)*dx,(Nz-1)*dx])

test.testT(t, 1)
test.testT(tc, 2)

#==========================================
# Initialisation du cas-test
#==========================================
Pref = 101320
Rhoref = 1.1765
gamma = 1.4
Rgp = 287.053
Tref = Pref/(Rhoref*Rgp)
c0 = numpy.zeros(1,dtype=numpy.float64)
c0[0] = numpy.sqrt(gamma*Rgp*Tref)

mach = 0.1
u0 = mach*c0[0]
v_adim=1./(numpy.sqrt(3.,dtype=numpy.float64)*c0)
v_adiminv = numpy.sqrt(3.,dtype=numpy.float64)*c0

# PARAMETRES DU VORTEX
R0 = L/10.
epsilon = 0.07*c0[0]
eps_2 = epsilon**2
dt = dx / (3. ** 0.5 * c0[0])
R0_2Inv= 1./(R0*R0)
R0_3Inv= 1./(R0*R0*R0)
x_0 = 0.9; y_0 = 0.5

C._initVars(t,'{centers:r2}=({centers:CoordinateX}-%20.16g)**2+({centers:CoordinateY}-%20.16g)**2'%(x_0,y_0))
C._initVars(t,'{centers:Sxx} = %20.16g*({centers:CoordinateY}-%20.16g)*%20.16g*({centers:CoordinateX}-%20.16g)*%20.16g*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,y_0,1./R0,x_0,R0_2Inv,R0_2Inv))
C._initVars(t,'{centers:Sxy} = 0.5*%20.16g*%20.16g*exp(-0.5*{centers:r2}*%20.16g)*(({centers:CoordinateY}-%20.16g)**2 - ({centers:CoordinateX}-%20.16g)**2)'%(epsilon,R0_3Inv,R0_2Inv,y_0,x_0))
C._initVars(t,'{centers:Sxz} = 0.0')
C._initVars(t,'{centers:Syy} = -%20.16g*({centers:CoordinateY}-%20.16g)*%20.16g*({centers:CoordinateX}-%20.16g)*%20.16g*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,y_0,R0_2Inv,x_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:Syz} = 0.0')
C._initVars(t,'{centers:Szz} = 0.0')
C._rmVars(t,['centers:r2'])

for i,z in enumerate(Internal.getZones(t)):
    if SOLVER_TYPE[i]==0:
       C._initVars(z,"{centers:r2}=({centers:CoordinateX}-%20.16g)**2+({centers:CoordinateY}-%20.16g)**2"%(x_0,y_0))
       C._initVars(z,'{centers:Density}= %20.16g*exp(-%20.16g/(2.*%20.16g**2)*exp(-{centers:r2}*%20.16g))'%(Rhoref,eps_2,c0[0],R0_2Inv))
       C._initVars(z,'{centers:VelocityX}=%20.16g-%20.16g*(({centers:CoordinateY}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(u0,epsilon,y_0,1./R0,R0_2Inv))
       C._initVars(z,'{centers:VelocityY}=%20.16g*(({centers:CoordinateX}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,x_0,1./R0,R0_2Inv))
       C._initVars(z,'{centers:VelocityZ}=0.')
       C._initVars(z,'{centers:Temperature}=1.')

       #ADIMENSIONNEMENT DE LA LBM
       C._initVars(z,'{centers:Density}={centers:Density}*%20.16g'%(1./Rhoref))
       C._initVars(z,'{centers:VelocityX}={centers:VelocityX}*%20.16g'%v_adim[0])
       C._initVars(z,'{centers:VelocityY}={centers:VelocityY}*%20.16g'%v_adim[0])
       C._initVars(z,'{centers:Sxx} = {centers:Sxx}*%20.16g'%dt)
       C._initVars(z,'{centers:Sxy} = {centers:Sxy}*%20.16g'%dt)
       C._initVars(z,'{centers:Syy} = {centers:Syy}*%20.16g'%dt)

    else:
        I._initWissocq(z, position=(x_0,y_0), Gamma=0.07, MInf=0.1, loc='centers')

C._addState(t, MInf=mach)
C._addState(t, adim='dim3',UInf=mach*c0,PInf=101320,RoInf=1.1765,LInf=1.)

VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']

#==========================================
# Preparation du calcul
#==========================================
ncyc = 0.5
#nit = int(ncyc*n/(1./numpy.sqrt(3.,dtype=numpy.float64)*mach))
nit = 100
print(nit)

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numb["omp_mode"]       = 0
numz = {}
numz["time_step"]          = dt
print(numz["time_step"])
numz["scheme"]             = "ausmpred"
numz["io_thread"]          = -1
# numz["LBM_HRR_sigma"]   = 0.
#numz['cache_blocking_I']=1000000
#numz['cache_blocking_J']=1000000
#numz['cache_blocking_K']=1000000

numz["LBM_taug"]        = 0.5#+nu_d/(c0[0]*c0[0]*numz["time_step"])
numz["LBM_collision"]   = "HRR"
numz["LBM_HRR_sigma"]   = 0.995
numz["LBM_NS"]          = 1

Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

#C.convertPyTree2File(t, 't.cgns')
#C.convertPyTree2File(tc, 'tc.cgns')
#test.testT(t,  7)
#test.testT(tc, 8)
#import sys;
#sys.exit()

print('avt warmup')
(t, tc, metrics) = Fast.warmup(t, tc, verbose=1)

'''
c= 0
for met in metrics:
  for me in met:
     print 'c=',c, me

  c+=1
'''
#import sys; sys.exit()
#C.convertPyTree2File(tc, 'tc.cgns')
#import sys; sys.exit()
#t2 = Internal.copyTree(t)
#test.testT(t2, 3)
#t2 = Internal.copyTree(tc)
#test.testT(t2, 4)
'''
t2 = Internal.copyTree(t)
test.testT(t2, 3)
import sys; sys.exit()
test.testT(tc, 4)
'''

#==========================================
# CALCUL
#==========================================
start = time.time()

tab_cons_NS = []
tab_cons_LBM = []

for it in range(1, nit+1):
   print('--------- iteration %d ---------'%it)
   #t2 = Internal.copyTree(t)
   #test.testT(t2, 6) 
   Fast._compute(t, metrics, it, tc, layer='Python')
   #t2 = Internal.copyTree(t)
   #test.testT(t2, 1000+it) 

#import sys; sys.exit()

#C.convertPyTree2File(t, 'tout.cgns')
test.testT(t, 5)
end = time.time()
print('')
print('Temps d''execution :',end-start)
print('')

'''plt.plot(tab_cons_NS+tab_cons_LBM)
plt.show()'''
#==========================================
# POST TRAITEMENT
#==========================================
tempVar = ['Q'+str(i) for i in range(1,20)]
C._rmVars(t, tempVar)
tempVar=['Q'+str(i)+'_M1' for i in range(1,20)]
C._rmVars(t,tempVar)

Internal.rmNodesByName(t, '.Solver#Param')
Internal.rmNodesByName(t, '.Solver#ownData')

tempVar = ['Density_P1','VelocityX_P1','VelocityY_P1','VelocityZ_P1','Temperature_P1']
C._rmVars(t, tempVar)
tempVar = ['Density_M1','VelocityX_M1','VelocityY_M1','VelocityZ_M1','Temperature_M1']
C._rmVars(t, tempVar)

for i,z in enumerate(Internal.getZones(t)):
   if SOLVER_TYPE[i]==0:
      C._initVars(z,'{centers:Density}={centers:Density}*%20.16g'%(Rhoref))
      C._initVars(z,'{centers:VelocityX} = {centers:VelocityX}*%20.16g'%(v_adiminv[0]))
      C._initVars(z,'{centers:VelocityY} = {centers:VelocityY}*%20.16g'%(v_adiminv[0]))
      C._initVars(z,'{centers:VelocityZ} = {centers:VelocityZ}*%20.16g'%(v_adiminv[0]))
      C._initVars(z,'{centers:Temperature} = {centers:Temperature}*%20.16g'%(Tref))
      C._initVars(z,'{centers:Sxx} = {centers:Sxx}*%20.16g'%(1./dt))
      C._initVars(z,'{centers:Sxy} = {centers:Sxy}*%20.16g'%(1./dt))
      C._initVars(z,'{centers:Syy} = {centers:Syy}*%20.16g'%(1./dt))
   i = i+1

for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO,storage=1)

# tempVar = ['Sxx','Sxy','Sxz','Syy','Syz','Szz']
# C._rmVars(t,tempVar)

VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO,storage=1)

t = P.computeGrad2(t,'centers:VelocityX')
t = P.computeGrad2(t,'centers:VelocityY')

VARSMACRO = ['gradxVelocityY','gradyVelocityX','gradxVelocityX','gradyVelocityY']
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO,storage=1)

C._initVars(t,'{centers:Vorticity} = {gradxVelocityY}-{gradyVelocityX}')
C._initVars(t,'{centers:DivV} = {gradxVelocityX}+{gradyVelocityY}')
tempVar=['gradxVelocityX','gradyVelocityX','gradzVelocityX']
C._rmVars(t,tempVar)
tempVar=['gradxVelocityY','gradyVelocityY','gradzVelocityY']
C._rmVars(t,tempVar)

VARSMACRO = ['Vorticity','DivV']
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO,storage=1)

Internal._rmGhostCells(t,t,2)
#C.convertPyTree2File(t, 'tout.cgns')
test.testT(t, 6)
