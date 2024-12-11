''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''                     CONVECTED VORTEX TEST CASE                           '''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import FastC.PyTree as FastC
import FastASLBM.PyTree as FastLBM
import Fast.PyTree as Fast
#import FastLBM.PyTree as FastLBM
import Connector.PyTree as X
import Transform.PyTree as T
import Post.PyTree as P
import KCore.test as test
import Apps.Fast.LBM as App

import math
import numpy
import time

myApp = App.LBM(format='single')

#==========================================
# Generation du maillage
#==========================================
x_min = 0.0
L = 1.
Nx=101; Ny = Nx; Nz = 11
#Nx=41; Ny = Nx; Nz = 3
n = Nx-1
dx = L/n
print(dx)
z_min = -4*dx
a1 = G.cart((x_min,x_min,z_min), (dx,dx,dx), (Nx,Ny,Nz))
t = C.newPyTree(['Base',a1])

t,tc = myApp.prepare(t, t_out=None, tc_out=None, NP=0, translation=[(Nx-1)*dx, (Ny-1)*dx,(Nz-1)*dx],NG=2)

test.testT(t,1)
test.testT(tc,2)
'''
'''

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
epsilon = 0.07*c0[0] #0.5*c0[0]
eps_2 = epsilon**2
dt = dx / (3. ** 0.5 * c0[0])
R0_2Inv= 1./(R0*R0)
R0_3Inv= 1./(R0*R0*R0)
x_0 = 0.5; y_0 = 0.5

# INITIALISATION
C._initVars(t,"{centers:r2}=({centers:CoordinateX}-%20.16g)**2+({centers:CoordinateY}-%20.16g)**2"%(x_0,y_0))
C._initVars(t,'{centers:Density}= %20.16g*exp(-%20.16g/(2.*%20.16g**2)*exp(-{centers:r2}*%20.16g))'%(Rhoref,eps_2,c0[0],R0_2Inv))
C._initVars(t,'{centers:VelocityX}=%20.16g-%20.16g*(({centers:CoordinateY}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(u0,epsilon,y_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:VelocityY}=%20.16g+%20.16g*(({centers:CoordinateX}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(0.,epsilon,x_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:VelocityZ}=0.')
C._initVars(t,'{centers:Temperature}=%20.16g'%(Tref))

C._initVars(t,'{centers:Sxx} = %20.16g*({centers:CoordinateY}-%20.16g)*%20.16g*({centers:CoordinateX}-%20.16g)*%20.16g*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,y_0,1./R0,x_0,R0_2Inv,R0_2Inv))
C._initVars(t,'{centers:Sxy} = 0.5*%20.16g*%20.16g*exp(-0.5*{centers:r2}*%20.16g)*(({centers:CoordinateY}-%20.16g)**2 - ({centers:CoordinateX}-%20.16g)**2)'%(epsilon,R0_3Inv,R0_2Inv,y_0,x_0))
C._initVars(t,'{centers:Syy} = -%20.16g*({centers:CoordinateY}-%20.16g)*%20.16g*({centers:CoordinateX}-%20.16g)*%20.16g*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,y_0,R0_2Inv,x_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:Sxz} = 0.0')
C._initVars(t,'{centers:Syz} = 0.0')
C._initVars(t,'{centers:Szz} = 0.0')

# ADIMENSIONNEMENT DE LA LBM
C._initVars(t,'{centers:Density}={centers:Density}*%20.16g'%(1./Rhoref))
C._initVars(t,'{centers:VelocityX}={centers:VelocityX}*%20.16g'%v_adim[0])
C._initVars(t,'{centers:VelocityY}={centers:VelocityY}*%20.16g'%v_adim[0])
C._initVars(t,'{centers:Sxx} = {centers:Sxx}*%20.16g'%dt)
C._initVars(t,'{centers:Sxy} = {centers:Sxy}*%20.16g'%dt)
C._initVars(t,'{centers:Syy} = {centers:Syy}*%20.16g'%dt)

C._rmVars(t,['centers:r2'])

C._addState(t, MInf=mach)
C._addState(t, adim='dim3', UInf=mach*c0, PInf=101320., RoInf=1.1765, LInf=1.)
VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature','Sxx','Sxy','Sxz','Syy','Syz','Szz']
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO)

#==========================================
# Preparation du calcul
#==========================================
ncyc = 0.01
nit = int(ncyc*n/(1./numpy.sqrt(3.,dtype=numpy.float64)*mach))
print(nit)

# Numerics
numb = {'temporal_scheme':'explicit', 'ss_iteration':20,'omp_mode':0}
numz = {'scheme':'ausmpred'}#,'io_thread':1}
numz["time_step"]=dt # pour l instant taug=0.5 (pas de viscosite)

# Helpful for debug
# numz['cache_blocking_I']=1000000
# numz['cache_blocking_J']=1000000
# numz['cache_blocking_K']=1000000

# LBM own parameters
nu_d = (15.6e-6)/Rhoref
numz["LBM_relax_time"]        = 0.5#+nu_d/(c0[0]*c0[0]*numz["time_step"])
numz["LBM_coll_model"]   = "HRR"
numz["LBM_hrr_sigma"]   = 0.98
# numz["LBM_HLBM"] = 0
# numz["LBM_dx"] = dx


myApp.set(numb=numb); myApp.set(numz=numz)
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

t = C.node2Center(t, ['CoordinateX','CoordinateY'])


#(t, tc, metrics)  = FastLBM.warmup(t, tc, flag_initprecise=1)
(t, tc, metrics)  = Fast.warmup(t, tc)

MESH = ['CoordinateX','CoordinateY']
for z in Internal.getZones(t):
    dim = Internal.getZoneDim(z)
    Nx = dim[1]-1-4; Ny = dim[2]-1-4; Nz = dim[3]-1-4
    GridCoordinates = numpy.zeros((Nx,Ny,2))
    gridcoord = Internal.getNodeFromName(z,'FlowSolution#Centers')
    for j,msh in enumerate(MESH):
        c = Internal.getNodeFromName(gridcoord,msh)[1]
        GridCoordinates[:,:,j] = c[2:-2,2:-2,0]

#test.testT(t,3)
#test.testT(tc,4)
#import os; sys.exit()


#==========================================
# CALCUL
#==========================================
tab_min = numpy.zeros((nit,2))
inst = 0.
start_all = time.time()

for it in range(0,nit+1):
    #inst = inst + dt
    print("--------- iteration %d -------------"%it)
    #start = time.time()
    #FastLBM._compute(t, metrics, it, tc, layer='c')
    #FastLBM._compute(t, metrics, it, tc, layer='Python')
    #Fast._compute(t, metrics, it, tc, layer='Python')
    Fast._compute(t, metrics, it, tc, layer='c')
    #end = time.time()
    #print('Temps d''execution :',end-start)
    # FastLBM.display_cpu_efficiency(t)

end_all = time.time()
print('')
print('Temps d''execution :',end_all-start_all)
print('')

#==========================================
# POST TRAITEMENT
#==========================================
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

#
VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature','Sxx','Sxy','Sxz','Syy','Syz','Szz']
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO)
#

#Internal._rmGhostCells(t,t,2)
#C.convertPyTree2File(t,"ASLBM.cgns")
test.testT(t,5)


#0.331663595885055
