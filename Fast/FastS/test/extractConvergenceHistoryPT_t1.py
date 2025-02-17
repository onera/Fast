# - extractConvergenceHistory (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import FastS.PyTree as FastS
import FastC.PyTree as FastC
import Converter.Internal as Internal
import Initiator.PyTree as I
import KCore.test as test
LOCAL = test.getLocal()

ni = 155; dx = 100./(ni-1); dz = 0.01
a1 = G.cart((-50,-50,0.), (dx,dx,dz), (ni,ni,2))
C._fillEmptyBCWith(a1, 'far', 'BCFarfield', dim=2)
I._initConst(a1, MInf=0.4, loc='centers')
C._addState(a1, 'GoverningEquations', 'Euler')
C._addState(a1, MInf=0.4)
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
       	FastS._displayTemporalCriteria(t, metrics, it, format='store')
    time += time_step

# Compute global residuals
FastS._calc_global_convergence(t)
# extraction des residus et creation du fichier "residus.dat"
FastS._extractConvergenceHistory(t, LOCAL+"/residus.dat")
test.testF(LOCAL+'/residus.dat', 1)
