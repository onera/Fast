## - Explicit local/local time stepping for flow past a square cylinder 3D -
## Note: mesh is "ugly"
import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree as Fast
import Fast.Utils as Utils
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import KCore.test as test
import Initiator.Adim as Adim
import Post.PyTree as P
import Transform.PyTree as T
import numpy
import sys
import os
import tarfile
import shutil

if not os.path.isfile('cgns_lts/t2DLTS.cgns'):
    tar = tarfile.open('cgns_lts.tar.gz', "r:gz")
    tar.extractall()
    tar.close()

t  = Fast.loadTree('cgns_lts/t2DLTS_ibm.cgns')
tc = Fast.loadTree('cgns_lts/tc2DLTS_ibm.cgns')
NIT                        = 100   # number of iterations
display_probe_freq         = 10    # iteration frequency to display modulo_verif

# Numerics
dt= 5.e-6
numb = {}
numb["temporal_scheme"]    = "explicit_local"
numb["modulo_verif"]       = display_probe_freq
numz = {}
numz["scheme"]             = "roe_min"
running_cfl                = 0.48
## Determines the maximum dt for the given flow field
dt_max = FastS.set_dt_lts(t,running_cfl,dt)
numz["time_step"]          = dt_max


it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)

# Numerics
FastC._setNum2Base(t, numb)
FastC._setNum2Zones(t, numz)
(t, tc, metrics) = FastS.warmup(t, tc) 

time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)
    
it_probe=0
for it in range(NIT):
    FastS._compute(t, metrics, it, tc,layer="Python")
    if it%display_probe_freq == 0:
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


                    

 


