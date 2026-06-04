# - Fast.IBM -
# Euler, para, frontType=1
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import FastC.PyTree as FastC
import FastS.Mpi as FastS
import Fast.IBM as App
import KCore.test as test

LOCAL = test.getLocal()

myApp = App.IBM(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3,
                "omp_mode":0})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

# Prepare
myApp.input_var.NP=Cmpi.size
t, tc = myApp.prepare('naca1DEuler.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns')
Internal._rmNodesFromName(tc,Internal.__GridCoordinates__)
Internal._rmNodesFromType(tc, 'Rind_t')
if Cmpi.rank == 0: test.testT(tc, 1)
Cmpi.barrier()

tc, graph = FastC.loadTree(LOCAL+'/tc.cgns', graph=True)

FastC._attributeNoPassTransfer(tc, graph=graph, cutoff=1.e-7, verbose=0)

Cmpi.convertPyTree2File(tc,LOCAL+'/tc.cgns')
tc, graph = FastC.loadTree(LOCAL+'/tc.cgns', graph=True)
t         = FastC.loadTree(LOCAL+'/t.cgns')

# Compute
numb={"temporal_scheme": "implicit", "ss_iteration":3, "omp_mode":0}
numz={"time_step": 0.0007, "scheme":"roe_min", "time_step_nature":"local", "cfl":4.}

it0 = 0.; time0 = 0.; NIT = 300
FastC._setNum2Base(t, numb); FastC._setNum2Zones(t, numz)

t, tc, metrics = FastS.warmup(t, tc, graph=graph, verbose=0)

it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)
time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)

if 'modulo_verif' in numb: moduloVerif = numb['modulo_verif']
else: moduloVerif = 200

for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph)
    if it%moduloVerif == 0:
        if Cmpi.rank == 0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it)
    time0 += time_step


Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)


FastC.save(t,LOCAL+'/restart.cgns')
if Cmpi.size > 1:
    Cmpi.barrier()

if Cmpi.rank == 0:
    t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
    Internal._rmNodesFromType(t, 'Rind_t')
    Internal._rmNodesFromName(t, '*M1')
    Internal._rmNodesFromName(t, '*P1')
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    test.testT(t, 2)
Cmpi.barrier()
