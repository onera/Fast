import Converter.PyTree as C
import Converter.Internal as Internal
import Apps.Fast.IBM as AppIBM
import Converter.Mpi as Cmpi
import Converter
import Fast.PyTree as Fast
import Generator.PyTree as G
import Connector.IBM as C_IBM
import Geom.PyTree as D
import Geom.IBM as G_IBM
import Transform.PyTree as T
import numpy
import KCore.test as test

import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS


STEP = 1
t_out = 't.cgns'
tc_out = 'tc.cgns'

snearsf = [0.035, 0.07]
frontType=1
#ext =2
#optimized=1
ext =1
optimized=-1
dfar=15
a = D.circle((0,0,0), 0.5, 0., 360.)
G_IBM._setDfar(a, dfar)
G_IBM._setIBCType(a, "Musker")
G_IBM._setSnear(a,0.01)

tb = C.newPyTree(["BODY", a])

#===============================================================================
# Etat de reference
#===============================================================================
Pref = 101320
Rhoref = 1.1765
gamma = 1.4
Rgp = 287.053
Tref = Pref/(Rhoref*Rgp)
c0 = numpy.zeros(1,dtype=numpy.float64)
c0[0] = numpy.sqrt(gamma*Rgp*Tref)

mach = 0.1
v_adim=1./(numpy.sqrt(3.,dtype=numpy.float64)*c0)
v_adiminv = numpy.sqrt(3.,dtype=numpy.float64)*c0

U0= mach*c0[0]
P0= Pref
L0= 1.
R0= Rhoref
T0= Tref
IFLOW='NSLaminar'

for base in Internal.getBases(tb):
  C._addState(base, 'GoverningEquations', IFLOW)
  C._addState(base, 'EquationDimension', 2)
  C._addState(base, UInf=U0, RoInf=R0, PInf=P0, LInf=L0, alphaZ=0., adim='dim3')

t,tc=AppIBM.prepare1(tb, t_out, tc_out,   vmin=15, snearsf=snearsf, frontType=frontType, cleanCellN=False,
                     nature=1, order=2, ext=ext, optimized=optimized, check=True)

C_IBM._buildConservativeFlux(t, tc, verbose=1)

t1c= Internal.copyRef(tc)
test.testT(t1c, 1)
t1= Internal.copyRef(t)
test.testT(t1, 2)

modulo_verif=10
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 30
numb["modulo_verif"]       = modulo_verif
numb["omp_mode"]           = 1
numz = {}
numz["scheme"]             = "ausmpred"
numz["slope"]             = "o3"
numz["ssdom_IJK"]          = [1000,1000,1000]
numz["time_step"]          = 0.0001

Fast._setNum2Base(t, numb)
Fast._setNum2Zones(t, numz)

refState = Internal.getNodeFromName(t, 'ReferenceState')
Mus = Internal.getNodeByName(refState, 'Mus')[1]
Mus[0] = (U0*Rhoref)/150                                                       # Re = 150

(t, tc, metrics) = FastS.warmup(t, tc, verbose=1 )

t1c= Internal.copyRef(tc)
test.testT(t1c, 3)
t1= Internal.copyRef(t)
test.testT(t1, 4)

for it in range(100):
     FastS._compute(t, metrics, it, tc)
     if it%(modulo_verif) == 0:
         FastS.displayTemporalCriteria(t, metrics, it)

test.testT(t, 5)

