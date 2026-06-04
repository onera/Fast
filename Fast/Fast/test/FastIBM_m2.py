# - Fast.IBM -
# Euler, para, frontType=1
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import FastC.PyTree as FastC
import FastS.Mpi as FastS
import Fast.IBM as App
import KCore.test as test
import Geom.PyTree as D

import Post.IBM as P_IBM

test.TOLERANCE = 1.e-8
LOCAL = test.getLocal()

FILEB = LOCAL+"/case.cgns"
FILEC = LOCAL+"/tc_restart.cgns"
FILED = LOCAL+"/tcw.cgns"

myApp = App.IBM(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":1,
                "omp_mode":0,
                "modulo_verif":200})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})
# case
a = D.sphere6((0.,0.,0.),1.,N=20)
for z in a:
    App._setSnear(z, 0.1)
    App._setDfar(z, 10.)
    App._setIBCType(z, 'Musker')

tb = C.newPyTree(['Body',a])
for base in Internal.getBases(tb):
    fl = Internal.newFlowEquationSet(parent=base)
    gov = Internal.newGoverningEquations(value='NSTurbulent', parent=fl)
eqdim = Internal.createNode('EquationDimension', '"int"', value=3, children=[])
turbmod = Internal.createNode('TurbulenceModel', 'TurbulenceModel_t', value='OneEquation_SpalartAllmaras', children=[])
for node in Internal.getNodesByName(tb,'FlowEquationSet'):
    Internal.addChild(node, eqdim)
    Internal.addChild(node, turbmod)

C._addState(tb, adim='adim1', MInf=0.1, alphaZ=0., alphaY=0., ReInf=40000., MutSMuInf=0.1, TurbLevelInf=1.e-4)

if Cmpi.rank == 0: C.convertPyTree2File(tb, FILEB)
Cmpi.barrier()

# Prepare
myApp.input_var.NP=Cmpi.size
myApp.input_var.vmin=11
t,tc = myApp.prepare(FILEB, t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns')
Internal._rmNodesFromType(tc, 'Rind_t')
Internal._rmNodesFromName(tc, Internal.__GridCoordinates__)
if Cmpi.rank == 0: test.testT(tc, 1)
Cmpi.barrier()

tc, graph = FastC.loadTree(LOCAL+'/tc.cgns', graph=True)

FastC._attributeNoPassTransfer(tc, graph=graph, cutoff=1.e-7, verbose=0)

Cmpi.convertPyTree2File(tc,LOCAL+'/tc.cgns')
tc, graph = FastC.loadTree(LOCAL+'/tc.cgns', graph=True)
t         = FastC.loadTree(LOCAL+'/t.cgns')

# case
# Compute
moduloVerif = 200
numb={"temporal_scheme": "implicit", "ss_iteration":1, "omp_mode":0, "modulo_verif":moduloVerif }
numz={"time_step": 0.0007, "scheme":"roe_min", "time_step_nature":"local", "cfl":4.}

FastC._setNum2Base(t, numb); FastC._setNum2Zones(t, numz)

t, tc, metrics = FastS.warmup(t, tc, graph=graph, verbose=0)

it0 = 0; time0 = 0.; NIT = 100
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)
time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)


for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph)
    if it%moduloVerif == 0:
        if Cmpi.rank == 0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it)
    time0 += time_step



# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)


FastC.save(t,LOCAL+'/restart.cgns')
if Cmpi.size > 1: Cmpi.barrier()

if Cmpi.rank == 0:
    t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
    Internal._rmNodesFromType(t, 'Rind_t')
    Internal._rmNodesFromName(t, '*M1')
    Internal._rmNodesFromName(t, '*P1')
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    test.testT(t, 2)

graphIBCDPost, ts = P_IBM.prepareSkinReconstruction(tb, tc, dimPb=3)
P_IBM._computeSkinVariables(ts, tc, graphIBCDPost, dimPb=3)
Cmpi.barrier()
if Cmpi.rank == 0: test.testT(ts, 3)
