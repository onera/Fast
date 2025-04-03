# - Fast.IBM -
# NS, para, frontType=2
import Fast.IBM as App
import Fast.FastIBM as FastIBM
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Converter.PyTree as C
import Converter.Internal as Internal
import math
import numpy as np
import KCore.test as test
test.TOLERANCE = 1.e-6

LOCAL = test.getLocal()

NEW = True

interpDataType = 1 # on suppose maillage non cartesion pour interp
cartesian = False

tb_2d = C.convertFile2PyTree('naca1DNS.cgns')

# Prepare
if NEW: t_2d, tc_2d = FastIBM.prepareIBMData('naca1DNS.cgns', None, None, check=False, frontType=2, cleanCellN=False)
else: t_2d, tc_2d = App.prepare1('naca1DNS.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', frontType=2, cleanCellN=False)

Internal._rmNodesFromType(tc_2d,'Rind_t')
Internal._rmNodesFromName(tc_2d,Internal.__GridCoordinates__)

####
# The following lines are to avoid regression
####
if NEW:
    floweq = Internal.getNodeByType(tb_2d, 'FlowEquationSet_t')
    refste = Internal.getNodeByType(tb_2d, 'ReferenceState_t')
    Internal._rmNodesByType(tc_2d, 'FlowEquationSet_t')
    Internal._rmNodesByType(tc_2d, 'ReferenceState_t')
    ibcdata = Internal.getNodeByName(tc_2d, '.Solver#IBCdefine')
    Internal._rmNodesByName(tc_2d, '.Solver#IBCdefine')
    base = Internal.getBases(tc_2d)[0]
    Internal.addChild(base, Internal.copyRef(floweq), -1)
    Internal.addChild(base, Internal.copyRef(refste), -1)
    Internal.addChild(base, ibcdata, -1)

    floweq = Internal.getNodeByType(tb_2d, 'FlowEquationSet_t')
    refste = Internal.getNodeByType(tb_2d, 'ReferenceState_t')
    Internal._rmNodesByType(t_2d, 'FlowEquationSet_t')
    Internal._rmNodesByType(t_2d, 'ReferenceState_t')
    base = Internal.getBases(t_2d)[0]
    Internal.addChild(base, Internal.copyRef(floweq), -1)
    Internal.addChild(base, Internal.copyRef(refste), -1)

    Internal._rmNodesByName(tc_2d, 'TurbulenceModel')
    Internal._rmNodesByName(t_2d, 'TurbulenceModel')

    for b in Internal.getBases(tc_2d):
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

test.testT(tc_2d, 1)
test.testT(t_2d, 2)

#on renomme les zones
i=0
for z in Internal.getZones(t_2d):
    z[0]= "Cart."+str(i)+"X0"
    i += 1

## determine dx=dy for each zone & store per zone
dict_ZEXT={}
hmin = 1.e30
for z in Internal.getZones(t_2d):
    h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
    print("dx=",h)
    dict_ZEXT[z[0]]=h
    if h < hmin : hmin = h

## go from dx to dx/dx_min
Nlevels=1
for i in dict_ZEXT:
    dict_ZEXT[i]= math.log( int(dict_ZEXT[i]/hmin + 0.00000001)  , 2)
    if dict_ZEXT[i] +1  > Nlevels : Nlevels = dict_ZEXT[i] +1

## get number of levels
print("Nlevel, hmin", Nlevels, hmin)

## Create the dict with Nz for each zone
dictNz      = {}

Nz_min = 8
Nz_max = 8
NzLoc=np.empty(int(Nlevels), np.int32)

for l in range( int(Nlevels) ):
    NzLoc[l] = Nz_max
Nlevels_tg = math.log( Nz_max/Nz_min, 2 ) +1

for z in Internal.getZones(t_2d):

    level = int( dict_ZEXT[ z[0] ] )
    print("Nz local",  z[0], level)
    dictNz[z[0]]=  NzLoc[level]
    print("Nz local", NzLoc[level], z[0], level)

extrusion ='cart'; span = 0.078
#extrusion ='cyl'; span = 22.5

if NEW: t_3d, tb_3d = FastIBM.extrudeCartesianZDir(t_2d, tb_2d, extrusion=extrusion, NPas=5, span=span, dictNz=dictNz, Ntranche=1)
else: t_3d, tb_3d = App.extrudeCartesian(t_2d, tb_2d, extrusion=extrusion, NPas=5, span=span, Ntranche=1, dictNz=dictNz)

for b in Internal.getBases(t_3d):
    if b[0] == 'CARTESIAN': b[0] = 'EXTRUSION' # Autrement Cartesian mode forced in Connector.Mpi

if NEW: t_3d, tc_3d = FastIBM.prepareIBMDataExtrude(tb_3d, None, None, t_3d, extrusion=extrusion, cartesian=cartesian, frontType=2)
else: t_3d, tc_3d = App.prepare1(tb_3d, None, None, t_in=t_3d, extrusion=extrusion, interpDataType=interpDataType, order=2, vmin=11, frontType=2)

####
# The following lines are to avoid regression
####
if NEW:
    floweq = Internal.getNodeByType(tb_2d, 'FlowEquationSet_t')
    refste = Internal.getNodeByType(tb_2d, 'ReferenceState_t')
    Internal._rmNodesByType(tc_3d, 'FlowEquationSet_t')
    Internal._rmNodesByType(tc_3d, 'ReferenceState_t')
    ibcdat = Internal.copyRef(Internal.getNodeFromName(tc_3d, '.Solver#IBCdefine'))
    Internal._rmNodesByName(tc_3d, '.Solver#IBCdefine')
    base = Internal.getBases(tc_3d)[0]
    Internal.addChild(base, Internal.copyRef(floweq), 0)
    Internal.addChild(base, Internal.copyRef(refste), 1)
    floweq = Internal.rmNodesByName(floweq, 'TurbulenceModel')
    Internal.addChild(base, Internal.copyRef(floweq), -1)
    Internal.addChild(base, Internal.copyRef(refste), -1)
    Internal.addChild(base, Internal.copyRef(ibcdat), -1)

    floweq = Internal.getNodeByType(tb_2d, 'FlowEquationSet_t')
    refste = Internal.getNodeByType(tb_2d, 'ReferenceState_t')
    Internal._rmNodesByType(t_3d, 'FlowEquationSet_t')
    Internal._rmNodesByType(t_3d, 'ReferenceState_t')
    ibcdat = Internal.copyRef(Internal.getNodeFromName(t_3d, '.Solver#IBCdefine'))
    Internal._rmNodesByName(t_3d, '.Solver#IBCdefine')
    base = Internal.getBases(t_3d)[0]
    Internal.addChild(base, Internal.copyRef(floweq), 0)
    Internal.addChild(base, Internal.copyRef(refste), 1)
    floweq = Internal.rmNodesByName(floweq, 'TurbulenceModel')
    Internal.addChild(base, Internal.copyRef(floweq), -1)
    Internal.addChild(base, Internal.copyRef(refste), -1)
    Internal.addChild(base, Internal.copyRef(ibcdat), -1)

    for b in Internal.getBases(tc_3d):
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

test.testT(tc_3d, 3)
test.testT(t_3d, 4)

t = Internal.copyRef(t_3d)
tc = Internal.copyRef(tc_3d)

# Compute
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 3
numb["omp_mode"]           = 1

numz = {}
numz["time_step"]          = 0.0007
numz["time_step_nature"]   = "local"
numz["cfl"]                = 4.
numz["scheme"]             = "roe_min"

it0 = 0.; time0 = 0.; NIT = 100
Fast._setNum2Base(t, numb); Fast._setNum2Zones(t, numz)

t, tc, metrics = FastS.warmup(t, tc)

time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)

for it in range(NIT):
    if it%100 == 0: print("it %d / %d"%(it, NIT), flush=True)
    FastS._compute(t, metrics, it, tc)
    time0 += time_step

Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
Fast.saveTree(t, LOCAL+'/restart.cgns', split='single', compress=0)
Fast.saveTree(tc, LOCAL+'/tc_restart.cgns', split='single')
t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(t, '.Solver#dtloc')
Internal._rmNodesFromType(t, 'Rind_t')
test.testT(t, 5)

# Post
