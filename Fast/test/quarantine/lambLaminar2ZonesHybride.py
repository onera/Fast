# - compute (pyTree) -
# Lamb vortex
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Converter.GhostCells as GC
import KCore.test as test
import Fast.PyTree as Fast
#import FastP.PyTree as FastP
import Initiator.PyTree as I
import Connector.PyTree as X
import Post.PyTree as P
import CPlot.PyTree as CPlot
import sys, numpy

t = C.convertFile2PyTree('tIvan.cgns')
tc= C.convertFile2PyTree('tcIvan.cgns')
mach = 0.7
#C.convertPyTree2File(t, 'tIvan.adf')
#C.convertPyTree2File(tc, 'tcIvan.adf')

#sys.exit()


C._addState(t, 'GoverningEquations', 'NSLaminar')
C._addState(t, MInf=mach)


I._initLamb(t, position=(45.,25.), Gamma=2., MInf=mach, loc='centers')



# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numb["modulo_verif"]       = 1
numz = {}
numz["time_step"]          = 0.03
numz["scheme"]             = "ausmpred"

Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

(t, tc, metrics) = Fast.warmup(t,tc)

zones = Internal.getZones(t)
for z in zones:
 sol = Internal.getNodeFromName1(z, 'FlowSolution#Centers')
 dens= Internal.getNodeFromName1(sol, 'Density')[1]
 dens[:]=dens[:]-0.10387

# Time Steps
nit = 300  ; time = 0.
for it in range(nit): 
    Fast._compute(t, metrics, it, tc)
    
    if it%1==0:
       Fast.display_temporal_criteria(t, metrics, it)
    
C.convertPyTree2File(t, 'out.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
