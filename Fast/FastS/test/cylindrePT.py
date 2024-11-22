# - compute (pyTree) -
# - Multibloc cylinder -
import Generator.PyTree as G
import Converter.PyTree as C
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Connector.PyTree as X
import Converter.Internal as Internal
import Transform.PyTree as T
import KCore.Adim as Adim

varsN = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']

MInf = 0.5
adim = Adim.adim1(MInf=MInf)

# Fichier classique multibloc
t = C.convertFile2PyTree('cylindreIn.cgns')
t = Internal.addGhostCells(t, t, 2, adaptBCs=1)
t = C.rmBCOfType(t, 'BCMatch')
t = C.fillEmptyBCWith(t, 'ov', 'BCOverlap', dim=2)
t = T.addkplane(t)

# Chimera
t = X.applyBCOverlaps(t, depth=2)
tc = C.node2Center(t)
tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse',
                     sameName=1, penalty=1, method='leastsquares', dim=2)
C._rmVars(tc, 'cellN') # tres important pour l'instant
#C.convertPyTree2File(tc, 'centers.cgns')
C._addState(t, MInf=MInf)
C._addState(t, 'GoverningEquations', 'Euler')
I._initConst(t, MInf=MInf, loc='centers')

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numz = {}
numz["time_step"]          = 0.003
numz["scheme"]             = "ausmpred"
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

(t, tc, metrics)  = FastS.warmup(t, tc)

nit = 3000; time = 0.
for it in range(nit):
    FastS._compute(t, metrics, it, tc)
    if it%20 == 0:
        print('- %d - %g'%(it, time))
        CPlot.display(t, dim=2, mode=3, scalarField=1)
    time += numz['time_step']

C.convertPyTree2File(t, 'out.cgns')
