# - compute (pyTree) -
# - Monobloc laminar flat plate -
import Converter.PyTree as C
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import FastC.PyTree as FastC
import FastS.PyTree as FastS

MInf = 0.2
t = C.convertFile2PyTree('flatPlate.cgns')
t = C.addState(t, 'GoverningEquations', 'NSLaminar')
t = C.addState(t, MInf=MInf)
t = I.initConst(t, MInf=MInf, loc='centers')
t = C.initVars(t, 'centers:ViscosityEddy', 0.)

numb = {'temporal_scheme':'explicit', 'ss_iteration':20}
numz = {'time_step':0.00002, 'scheme':'ausmpred'}
FastC._setNum2Zones(t, numz) ; FastC._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

nit = 30000 ; time = 0.

for it in range(nit):
    FastS._compute(t, metrics, it)
    if it%20 == 0:
        print('- %d - %g'%(it, time))
        CPlot.display(t, dim=2, mode=3, scalarField=2)
    time += numz['time_step']

C.convertPyTree2File(t, 'out.cgns')
