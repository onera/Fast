# - convergenceHistory (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Connector.PyTree as X
import Converter.Internal as Internal
import Initiator.PyTree as I

ni = 155; dx = 100./(ni-1); dz = 0.01
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
a1 = C.fillEmptyBCWith(a1, 'far', 'BCFarfield', dim=2)
a1 = I.initConst(a1, MInf=0.4, loc='centers')
a1 = C.addState(a1, 'GoverningEquations', 'Euler')
a1 = C.addState(a1, MInf=0.4)
t = C.newPyTree(['Base', a1])

# Numerics
modulo_verif = 20
numb = {}
numb["temporal_scheme"]    = "implicit"
numb["ss_iteration"]       = 3
numb["modulo_verif"]       = modulo_verif
numz = {}
numz["time_step"]          = 0.0007
numz["time_step_nature"]   = "local"
numz["cfl"]                = 4.0
numz["scheme"]             = "roe"
numz["slope"]              = "minmod"
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, None)

# Number or records to store residuals 
nrec = 100//modulo_verif

#To remove old ConvergenceHistory nodes 
C._rmNodes(t, "ZoneConvergenceHistory")

#Convergence history with nrec records
FastS._createConvergenceHistory(t, nrec)

nit = 100; time = 0
time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)
for it in range(nit):
    FastS._compute(t, metrics, it, tc)
    if it%modulo_verif == 0:
        FastS.displayTemporalCriteria(t, metrics, it, format='store')
    time += time_step

# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=nit)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time)
C.convertPyTree2File(t, 'out.cgns')
