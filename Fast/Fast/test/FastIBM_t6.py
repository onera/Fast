# - FastIBM -
# NS, para, frontType=1
import Apps.Fast.FastIBM as FastIBM
import FastS.PyTree as FastS
import Fast.PyTree as Fast
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
test.TOLERANCE = 1.e-7

LOCAL = test.getLocal()

# Prepare
t,tc = FastIBM.prepareIBMData('naca1DNS.cgns', None, None, vmin=21, frontType=1)
Internal._rmNodesFromType(tc, 'Rind_t')
Internal._rmNodesFromName(tc, Internal.__GridCoordinates__)
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
Internal._rmNodesByName(t, 'TurbulenceModel')

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

test.testT(tc, 1)

# Compute
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 3
numb["omp_mode"]           = 0

numz = {}
numz["time_step"]          = 0.0007
numz["time_step_nature"]   = "local"
numz["cfl"]                = 4.
numz["scheme"]             = "roe_min"

it0 = 0.; time0 = 0.; NIT = 300
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
test.testT(t, 2)

# Post
