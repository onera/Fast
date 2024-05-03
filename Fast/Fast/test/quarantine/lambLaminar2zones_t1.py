# - compute (pyTree) -
# Lamb vortex
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Converter.GhostCells as GC
import KCore.test as test

import FastP.PyTree as FastP
import Fast.PyTree as Fast
import Initiator.PyTree as I
import Connector.PyTree as X
import Post.PyTree as P
import sys

a = G.cartNGon((0,0,0), (0.1,0.1,1), (180,180,3))
t = C.newPyTree(['Base',a])
C._fillEmptyBCWith(t, 'extrap', 'BCExtrapolate', dim=3)

#C.convertPyTree2File(t, 'out_avtsplit.cgns'
t = T.splitNParts(t, 2)

t = X.connectMatch(t,dim=3)

Internal._adaptNFace2PE(t, remove=False) 

#
ncouche =2
t = GC.addGhostCellsNG(t, nlayers=ncouche) 

mach = 0.7
C._addState(t, 'GoverningEquations', 'NSLaminar')
C._addState(t, MInf=mach, ReInf=30.)

tc = C.node2Center(t)
X._setInterpData(t,tc,loc='centers',storage='inverse',itype='abutting')

I._initLamb(t, position=(4.,4.), Gamma=2., MInf=mach, loc='centers')


# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["ss_iteration"]       = 20
numb["modulo_verif"]       = 5
numz = {}
numz["time_step"]         = 0.03
numz["scheme"]            = "ausmpred"
Fast._setNum2Zones(t, numz); Fast._setNum2Base(t, numb)

zones = Internal.getZones(t)
znew=[]
for z in zones:
   nfaceElt = Internal.getNodeFromName1(z, 'NFaceElements')
   Elt0     = Internal.getNodeFromName1(nfaceElt, 'IntExt')[1][0]
   Elt1     = Internal.getNodeFromName1(nfaceElt, 'IntExt')[1][1]
   Elt2     = Internal.getNodeFromName1(nfaceElt, 'IntExt')[1][2]
   Elt3     = Internal.getNodeFromName1(nfaceElt, 'IntExt')[1][3]
   Elt4     = Internal.getNodeFromName1(nfaceElt, 'IntExt')[1][4]
   print('nb elt 0  1  2 3 4= ', Elt0, Elt1, Elt2, Elt3,Elt4, z[0])
   Index = list(range(Elt0+Elt1)) + list(range(Elt0+Elt1+Elt3,Elt0+Elt1+Elt3+ Elt2 ))
   znew.append( T.subzone(z, Index , type='elements'))


#(t, tc, metrics) = FastP.warmup(t, tc)
(t, tc, metrics) = Fast.warmup(t, tc)


zones = Internal.getZones(t)


param_int  = Internal.getNodeFromName2(zones[0], 'Parameter_int')[1]

C._initVars(t, '{centers:VelocityY} = {centers:VelocityY} + 1.')

#recompactage obligatoire, car _initvar is not in place

# Time Steps
nit =  25  ; time = 0.
for it in range(nit): 
    Fast._compute(t, metrics,  nit, tc,layer='C')

    if it%5==0:
       Fast.display_temporal_criteria(t, metrics, it)

C.convertPyTree2File(t, 'out.cgns')
#sys.exit()

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)

