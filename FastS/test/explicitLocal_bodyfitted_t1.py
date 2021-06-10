## - Explicit local/local time stepping for flow past a naca 0012 2D -
## Note: mesh is "ugly"
import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree  as Fast
import Fast.Utils as Utils
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import KCore.test as test
import Transform.PyTree as T
import math
import numpy
import sys
import os
import tarfile
import shutil

## BEGINNING OF COMPUTE
if not os.path.isfile('cgns_lts/t3DLTS.cgns'):
    tar = tarfile.open('cgns_lts.tar.gz', "r:gz")
    tar.extractall()
    tar.close()

t  = Fast.loadTree('cgns_lts/t3DLTS.cgns')
tc = Fast.loadTree('cgns_lts/tc3DLTS.cgns')

NIT                        = 100     # number of iterations
display_probe_freq         = 10      # iteration frequency to display modulo_verif


## Setting dt for simulation
numb = {}
numb["temporal_scheme"]    = "explicit_local"
numb["modulo_verif"]       = display_probe_freq
numz = {}
numz["scheme"]             = "ausmpred"
running_cfl                = 0.25
dt_max = FastS.set_dt_lts(t,running_cfl)
numz["time_step"]          = dt_max

it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)

# Numerics
FastC._setNum2Base(t, numb)
## Create for CGNSBase_t
## .Solver#define

FastC._setNum2Zones(t, numz)
## Create for each Zone_t
## .Solver#define

print('++++ pre-warmup ++++')
(t, tc, metrics) = FastS.warmup(t, tc) 
print('++++ post-warmup ++++')

time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)
for it in range(NIT):
    FastS._compute(t, metrics, it, tc,layer='Python')   
    if it%display_probe_freq == 0:
        print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it, format='double')
        
    time0 += time_step

# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
shutil.rmtree("cgns_lts")
#C.convertPyTree2File(t, "out.cgns")



 


