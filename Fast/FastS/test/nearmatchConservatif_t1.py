import Converter.PyTree as C
import Converter.Internal as Internal
import Fast.PyTree as Fast
import Connector.IBM as X_IBM
import Geom.PyTree as D
import Geom.IBM as G_IBM
import numpy
import KCore.test as test

import Fast.PyTree as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS

STEP = 1

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

# t,tc=AppIBM.prepare1(tb, None, None,   vmin=15, snearsf=snearsf, frontType=frontType, cleanCellN=False,
#                      nature=1, order=2, ext=ext, optimized=optimized, check=False)

t,tc=X_IBM.prepareIBMData(tb, None, None, vmin=15, snearsf=snearsf, frontType=frontType,
                    ext=ext+1, optimized=optimized, check=False, cleanCellN=False)

####
# The following lines are to avoid regression
####
floweq = Internal.getNodeByType(tc, 'FlowEquationSet_t')
Internal._rmNodesByType(tc, 'FlowEquationSet_t')
refste = Internal.getNodeByType(tc, 'ReferenceState_t')
Internal._rmNodesByType(tc, 'ReferenceState_t')
ibcdata = Internal.getNodeByName(tc, '.Solver#IBCdefine')
Internal._rmNodesByName(tc, '.Solver#IBCdefine')
base = Internal.getBases(tc)[0]
Internal.addChild(base, floweq, -1)
Internal.addChild(base, refste, -1)
Internal.addChild(base, ibcdata, -1)

Internal.getNodeFromName(tc, 'CARTESIAN')[1][0] = 3
Internal._rmNodesByName(tc, 'TurbulenceModel')

floweq = Internal.getNodeByType(t, 'FlowEquationSet_t')
Internal._rmNodesByType(t, 'FlowEquationSet_t')
refste = Internal.getNodeByType(t, 'ReferenceState_t')
Internal._rmNodesByType(t, 'ReferenceState_t')
base = Internal.getBases(t)[0]
Internal.addChild(base, floweq, -1)
Internal.addChild(base, refste, -1)

for z in Internal.getZones(t):
    TurbulentDistance = Internal.getNodeFromName(z, 'TurbulentDistance')[1]
    if numpy.max(TurbulentDistance) > 1.2:
        Internal.getNodeFromName(z, 'TurbulentDistance')[1] = numpy.zeros_like(TurbulentDistance)
        Internal.getNodeFromName(z, 'TurbulentDistanceAllBC')[1] = numpy.zeros_like(TurbulentDistance)

for b in Internal.getBases(tc):
    for z in Internal.getZones(b):
        pos = 0
        z2 = Internal.copyRef(z)
        dictOfIBCDZones = {zs[0]: Internal.copyRef(zs) for zs in z2[2] if 'IBCD' in zs[0]}
        dictOfIBCRZones = {zs[0].split('_')[-1]: [] for zs in z2[2] if 'IBCD' in zs[0]}
        rcvzone = ''
        for zs in z2[2]:
            if 'IBCD' in zs[0]:
                rcvzoneL = zs[0].split('_')[-1]
                dictOfIBCRZones[rcvzoneL].append(zs[0])
                listOfIBCRZonesL = list(dictOfIBCRZones[rcvzoneL])
                pos += 1
                for zname in dictOfIBCRZones[rcvzoneL]:
                    listOfIBCRZonesL.append(zname)
                    Internal.addChild(z, dictOfIBCDZones[zname], pos)
                    pos +=1
                dictOfIBCRZones[rcvzoneL] = list(listOfIBCRZonesL)
            else:
                pos += 1
####
                
X_IBM._buildConservativeFlux(t, tc, verbose=1)

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
