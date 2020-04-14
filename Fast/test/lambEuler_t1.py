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
import FastP.PyTree as FastP
import Initiator.PyTree as I
import Connector.PyTree as X
import Post.PyTree as P
import CPlot.PyTree as CPlot
import sys

a = G.cartNGon((0,0,0), (1,1,1), (100,100,3))
t = C.newPyTree(['Base',a])
C._fillEmptyBCWith(t, 'extrap', 'BCExtrapolate', dim=3)
Internal._adaptNFace2PE(t, remove=False) 

#
ncouche =2
t= GC.addGhostCellsNG(t, nlayers=ncouche) 

mach = 0.7

C._addState(t, 'GoverningEquations', 'Euler')
#C._addState(t, 'GoverningEquations', 'NSLaminar')
C._addState(t, MInf=mach)

I._initLamb(t, position=(25.,25.), Gamma=2., MInf=mach, loc='centers')

#tc = C.node2Center(t)
#tc = X.setInterpData(t, tc, nature=1, loc='centers', storage='inverse', 
#                     sameName=1, dim=3)
#C.convertPyTree2File(tc, 'tc.cgns')

# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
#numb["omp_mode"]          = 1
numz = {}
numz["time_step"]          = 0.3
numz["scheme"]             = "ausmpred"

Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

print('NEQQQQ=') 
#P._computeVariables(t, ['centers:VelocityX'  ])
(t, tc, metrics) = Fast.warmup(t)
#C.convertPyTree2File(t, 'outwithGH.cgns')
#sys.exit()

zones      = Internal.getNodesFromType2(t, 'Zone_t')
param_int  = Internal.getNodeFromName2(zones[0], 'Parameter_int')[1]
print('NEQQQQ=', param_int[36]) 

print("ti", metrics[0][1][1,0])
import numpy
ti = metrics[0][0]
vol= metrics[0][1]
size_ti   = numpy.shape( ti)
size_vol = numpy.shape(vol)
print("ti", size_ti)
print("vol", size_vol)
#for i in range(size_ti[0] ):
 #print "txyz", metrics[0][0][i,0],  metrics[0][0][i,1] ,metrics[0][0][i,2], i
 #print "fxyz", metrics[0][0][i,3],  metrics[0][0][i,4] ,metrics[0][0][i,5], i
#for i in range(size_vol[0] ):
 #print "vol", vol[i,0], vol[i,1],i
 #print "cxyz", vol[i,1],  vol[i,2] ,vol[i,3], i

#C.convertPyTree2File(t, 'out.cgns')
#sys.exit()


#C._initVars(t, '{centers:Density} = 1 + 0.0001 * {centers:CoordinateX}* {centers:CoordinateX}')
#C._initVars(t, '{centers:Temperature} = 1')
#C._initVars(t, '{centers:VelocityX} = 1')
#C._initVars(t, '{centers:VelocityY} = 1')
C._initVars(t, '{centers:VelocityY} = {centers:VelocityY} + 0.7')
#C.convertPyTree2File(t, "in.cgns")
#sys.exit()
# Time Steps
#nit =  50  ; time = 0.
nit =  50  ; time = 0.
for it in range(nit): 
    print('nit=',it)
    Fast._compute(t, metrics, nit)

#C.convertPyTree2File(t, 'out.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
