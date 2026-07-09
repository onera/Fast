import Fast.IBMO as App
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.test as test
test.TOLERANCE = 1.e-8

LOCAL = test.getLocal()

NP = Cmpi.size
rank = Cmpi.rank
NIT = 100
FILE = 'naca_IBMO.cgns'

myApp = App.IBMO(format='single')
myApp.set(numb={"temporal_scheme": "explicit",
                "ss_iteration":5,
                "omp_mode":0})
myApp.set(numz={"time_step": 0.002,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":0.5})

t,tc = myApp.prepare(FILE, t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', expand=3, vmin=11, check=False, NP=Cmpi.size, distrib=True)

cartBase = Internal.getNodeFromName(t,'CARTESIAN')
Internal._rmNodesFromType(cartBase,'Rind_t')
test.testT(t,1)

tc = C.convertFile2PyTree(LOCAL+'/tc.cgns')

FastC._attributeNoPassTransfer(tc, verbose=0, cutoff=1.e-12)



time_step= 0.002
numb={"temporal_scheme": "explicit", "ss_iteration":5, "omp_mode":0}
numz={"time_step": time_step, "scheme":"roe_min", "time_step_nature":"local", "cfl":0.5}
FastC._setNum2Base(t, numb); FastC._setNum2Zones(t, numz)

Internal._rmNodesByName(t,'senseurType')

t, tc, metrics = FastS.warmup(t, tc)

time0 = 0.
for it in range(NIT):
    FastS._compute(t, metrics, it, tc)
    time0 += time_step

Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
cartBase = Internal.getNodeFromName(t,'CARTESIAN')
Internal._rmNodesFromType(cartBase,'Rind_t')

####
# The following lines are to avoid regression since the removal of sortByName in FastS warmup
####
Internal._sortByName(t, recursive=False)
cgnslibver = Internal.getNodeByType(t, 'CGNSLibraryVersion_t')
Internal._rmNodesByType(t, 'CGNSLibraryVersion_t')
Internal.addChild(t, cgnslibver, 0)
####
test.testT(t,2)
