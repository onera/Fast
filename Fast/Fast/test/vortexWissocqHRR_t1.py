''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''                     CONVECTED VORTEX TEST CASE                           '''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import FastC.PyTree as FastC
import Fast.PyTree as Fast
import FastASLBM.PyTree as FastLBM
import Connector.PyTree as X
import Apps.Fast.LBM as Apps_LBM
import KCore.test as test

import numpy
import time

#==========================================
# Generation du maillage
#==========================================
x_min = 0.0; L = 1.
Nx=101; Ny = Nx; Nz = 11
n = Nx-1; dx = L/n; print('Grid spacing=',dx)
z_min = -4*dx

a1 = G.cart((x_min,x_min,z_min), (dx,dx,dx), (Nx,Ny,Nz))
t = C.newPyTree(['Base',a1])
t,tc = Apps_LBM.prepare(t, t_out=None, tc_out=None, NP=0, translation=[(Nx-1)*dx, (Ny-1)*dx,(Nz-1)*dx],NG=2)

test.testT(t,1); test.testT(tc,2)

#==========================================
# Initialisation du cas-test
#==========================================
Pref = 101320.
Rhoref = 1.1765
gamma = 1.4
Rgp = 287.053
Tref = Pref/(Rhoref*Rgp)
c0 = numpy.zeros(1,dtype=numpy.float64)
c0[0] = numpy.sqrt(gamma*Rgp*Tref)
mach = 0.1
u0 = mach*c0[0]
visco_nu = 0.

# Parametres du vortex
R0 = L/10.
epsilon = 0.07*c0[0]; eps_2 = epsilon**2
R0_2Inv= 1./(R0*R0); R0_3Inv= 1./(R0*R0*R0)
x_0 = 0.5; y_0 = 0.5

# Initialisation : champs macros
C._initVars(t,"{centers:r2}=({centers:CoordinateX}-%20.16g)**2+({centers:CoordinateY}-%20.16g)**2"%(x_0,y_0))
C._initVars(t,'{centers:Density}= %20.16g*exp(-%20.16g/(2.*%20.16g**2)*exp(-{centers:r2}*%20.16g))'%(Rhoref,eps_2,c0[0],R0_2Inv))
C._initVars(t,'{centers:VelocityX}=%20.16g-%20.16g*(({centers:CoordinateY}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(u0,epsilon,y_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:VelocityY}=%20.16g+%20.16g*(({centers:CoordinateX}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(0.,epsilon,x_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:VelocityZ}=0.')
C._initVars(t,'{centers:Temperature}=%20.16g'%(Tref))

# Initialisation : gradients
C._initVars(t,'{centers:Sxx} = %20.16g*({centers:CoordinateY}-%20.16g)*%20.16g*({centers:CoordinateX}-%20.16g)*%20.16g*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,y_0,1./R0,x_0,R0_2Inv,R0_2Inv))
C._initVars(t,'{centers:Sxy} = 0.5*%20.16g*%20.16g*exp(-0.5*{centers:r2}*%20.16g)*(({centers:CoordinateY}-%20.16g)**2 - ({centers:CoordinateX}-%20.16g)**2)'%(epsilon,R0_3Inv,R0_2Inv,y_0,x_0))
C._initVars(t,'{centers:Syy} = -%20.16g*({centers:CoordinateY}-%20.16g)*%20.16g*({centers:CoordinateX}-%20.16g)*%20.16g*exp(-0.5*{centers:r2}*%20.16g)'%(epsilon,y_0,R0_2Inv,x_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:Sxz} = 0.0')
C._initVars(t,'{centers:Syz} = 0.0')
C._initVars(t,'{centers:Szz} = 0.0')
C._rmVars(t,['centers:r2'])

C._addState(t, adim='dim3',UInf=u0,PInf=Pref,RoInf=Rhoref,LInf=L); C._addState(t, 'Mus', visco_nu)
dt_list, tau_list = FastLBM.getTimeStepAndRelaxTime(t); print(dt_list,tau_list)
# FastLBM._convertToLatticeUnits(t,dt_list)

#==========================================
# Parametres de calcul
#==========================================

# Parametres globaux
numb = {'temporal_scheme':'explicit', 'ss_iteration':20,'omp_mode':0}
numz = {'scheme':'ausmpred'}#,'io_thread':1}
numz["time_step"]=dt_list[0][0] # pour l instant taug=0.5 (pas de viscosite)
# Utile pour debug
# numz['cache_blocking_I']=1000000
# numz['cache_blocking_J']=1000000
# numz['cache_blocking_K']=1000000

# Parametres propres a la LBM
# numz["LBM_velocity_set"] = "D3Q19"
numz["LBM_coll_model"]   = "HRR"
numz["LBM_hrr_sigma"]   = 0.98
numz["LBM_relax_time"]   = tau_list[0][0]

FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

#==========================================
# Warmup
#==========================================
t = C.node2Center(t, ['CoordinateX','CoordinateY'])
(t, tc, metrics)  = Fast.warmup(t, tc)

#==========================================
# Calcul
#==========================================
ncyc = 0.01
nit = int(ncyc*n/(1./numpy.sqrt(3.,dtype=numpy.float64)*mach))
print(nit)

start_all = time.time()

for it in range(0,nit+1):
    print("--------- iteration %d -------------"%it)
    Fast._compute(t, metrics, it, tc)

end_all = time.time()
print('')
print('Temps d''execution :',end_all-start_all)
print('')

#==========================================
# POST TRAITEMENT
#==========================================

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(t, "*_P1")
Internal._rmNodesByName(t, "niveaux_temps")

FastLBM._convertToLatticeUnits(t,dt_list)

VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature','Sxx','Sxy','Sxz','Syy','Syz','Szz']
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO)

#Internal._rmGhostCells(t,t,2)
# C.convertPyTree2File(t,"ASLBM.cgns")
Internal._rmNodesByName(t, "*_P1")
test.testT(t,5)
