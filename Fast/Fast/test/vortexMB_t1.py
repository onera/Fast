# - PSEUDO ISENTROPIC VORTEX MULTIBLOCK (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Fast.PyTree as Fast
import FastLBM.PyTree as FastLBM
import Connector.PyTree as X
import Transform.PyTree as T
import KCore.test as test
import math
import numpy
import Apps.Fast.LBM as App

myApp = App.LBM(format='single')
nit = 300
VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']

# Geometry and mesh
x_min = -0.5
dx = 0.01
L = 1.
NG=1
Nx=101; Ny = Nx; Nz = 9
c0=343.2
n = Nx-1
tree=dict()
tree['mach']=0.1
tree['maille']=0.01
tree['timestep']=tree['maille']/(c0*numpy.sqrt(3.))
tree['reynolds']=100.
tree['char_length']=L/10. # R0
tree['rho0']=1.

# Fluid
rho0 = 1.
u0 = tree["mach"]* c0

# Vortex
R0 = tree['char_length']
kappa = 0.14
dt = dx / (3. ** 0.5 * c0)
#-------------------------------------------------------
#-------------------------------------------------------
z_min = -4*dx
a1 = G.cart((x_min,x_min,z_min), (dx,dx,dx), (Nx,Ny,Nz))
x_0 = C.getMeanValue(a1,'CoordinateX')
y_0 = C.getMeanValue(a1,'CoordinateY')
z_0 = C.getMeanValue(a1,'CoordinateZ')
zmean = C.getMeanValue(a1,'CoordinateZ')
a1 = T.splitNParts(a1,8)
t = C.newPyTree(['Base',a1])
t,tc = myApp.prepare(t, t_out=None, tc_out=None, NP=0, translation=[(Nx-1)*dx, (Ny-1)*dx,(Nz-1)*dx])

#-------------------------
# Initialization
#-------------------------
rho0Kap2 = rho0*kappa**2
kapC0=kappa*c0
R0_2Inv= 1./(R0*R0)
v_adim=1./(math.sqrt(3)*c0)
#
C._initVars(t,"{centers:r2}=({centers:CoordinateX}-%g)**2+({centers:CoordinateY}-%g)**2"%(x_0,y_0))
C._initVars(t,'{centers:Density}=1.-0.5*%g*exp(1.-{centers:r2}*%g)'%(kappa*kappa,R0_2Inv))
C._initVars(t,'{centers:VelocityX}=%g-%g*({centers:CoordinateY}*%g)*exp(0.5*(1.-{centers:r2}*%g))'%(u0,kappa*c0,1/R0,R0_2Inv))
C._initVars(t,'{centers:VelocityY}=%g*({centers:CoordinateX}*%g)*exp(0.5*(1.-{centers:r2}*%g))'%(kappa*c0,1/R0,R0_2Inv))
C._rmVars(t,["centers:r2"])
C._initVars(t,'{centers:Temperature}=1.')# pour l instant ne sert qu a passer dans le warmup, car isotherme
C._initVars(t,'{centers:VelocityZ}=0.')
C._initVars(t,"{centers:TurbulentDistance} =1.")
C._initVars(t,"{centers:cellN_IBC_LBM_1} =0.")
C._initVars(t,"{centers:cellN_IBC_LBM_2} =0.")
C._initVars(t,"{centers:cellN_IBC_LBM_3} =0.")
#-------------------------
# Adimensionnement
#-------------------------
C._initVars(t,'{centers:VelocityX}={centers:VelocityX}*%g'%v_adim)
C._initVars(t,'{centers:VelocityY}={centers:VelocityY}*%g'%v_adim)

for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO)
#-------------------------
# Compute
#-------------------------
# Numerics
numb = {'temporal_scheme':'explicit', 'ss_iteration':20,'omp_mode':0}
numz = {'scheme':'ausmpred'}
numz['cache_blocking_I']=1000000
numz['cache_blocking_J']=1000000
numz['cache_blocking_K']=1000000
numz['cache_blocking_I']=1000000
numz["time_step"]=dt # pour l instant pas de viscosite
myApp.set(numb=numb)
myApp.set(numz=numz)
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)
#C.convertPyTree2File(t,"restart.cgns")
(t, tc, metrics)  = FastLBM.warmup(t, tc)

for it in range(1,nit+1):
    print("--------- iteration %d -------------"%it)
    FastLBM._compute(t, metrics, it, tc,layer='Python',nittotal=nit)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')

# POST
Internal._rmNodesByName(t,'*M1*')
Internal._rmNodesByName(t,'*Sx*')
Internal._rmNodesByName(t,'*Sy*')
Internal._rmNodesByName(t,'*Sz*')
Internal._rmNodesByName(t,'*cell*')
Internal._rmNodesByName(t,'*Qstar*')
Internal._rmNodesByName(t,'*Qeq*')
Internal._rmNodesByName(t,'*Qneq*')
Internal._rmNodesByName(t,'*SpongeCoef')

# reconstruction des variables macro
#FastLBM.convert_to_dimensional(t,rho0,dx,dt)
v_adiminv=1./v_adim
C._initVars(t,"{centers:VelocityX}={centers:VelocityX}*%20.16g/{centers:Density}"%(v_adiminv))
C._initVars(t,'{centers:VelocityY}={centers:VelocityY}*%20.16g/{centers:Density}'%(v_adiminv))
C._initVars(t,'{centers:VelocityZ}={centers:VelocityZ}*%20.16g/{centers:Density}'%(v_adiminv))

# test.testT(t,1)
#t=FastLBM.rmGhostCells(t,NG)
#C.convertPyTree2File(t,"restart.cgns")
