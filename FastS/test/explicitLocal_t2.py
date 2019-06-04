# - compute (pyTree) -
# - Lamb vortex sur 2 domaines -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter as Co
import Generator.PyTree as G
import Initiator.PyTree as I
import CPlot.PyTree as CPlot
import Fast.PyTree as Fast_
import Fast.Internal as internal
import FastS.PyTree as FastS
import Post.PyTree as P
import Fast.Utils
import Connector.PyTree as X
import Transform.PyTree as T
import sys
import Converter.Internal as Internal
import numpy
import math
import timeit
import pickle
import time as Time
import KCore.test as test


#### Maillage ####
a1 = G.cart((0,0,0),(0.05,0.05,0.1), (351,201,1))
a1 = C.addVars(a1, 'centers:niveaux_temps')
a1 = C.initVars(a1, 'centers:niveaux_temps',1)

a2 = G.cart((17.5,0,0),(0.05,0.05,0.1), (151,201,1))
a2 = C.addVars(a2, 'centers:niveaux_temps')
a2 = C.initVars(a2, 'centers:niveaux_temps',2)

a3 = G.cart((0,10,0),(0.05,0.05,0.1), (351,201,1))
a3 = C.addVars(a3, 'centers:niveaux_temps')
a3 = C.initVars(a3, 'centers:niveaux_temps',2)

a4 = G.cart((17.5,10,0),(0.05,0.05,0.1), (151,201,1))
a4 = C.addVars(a4, 'centers:niveaux_temps')
a4 = C.initVars(a4, 'centers:niveaux_temps',4)

### Arbre t ####
t = C.newPyTree(['Base']) ; t[2][1][2] += [a1,a2,a3,a4]

### Conditions initiales ####
t = C.addState(t, 'GoverningEquations', 'Euler')
t = C.addState(t, 'EquationDimension', 2)
t = C.addState(t, MInf=0.4, ReInf=100)
t = I.initLamb(t, position=(15,10.), Gamma=2., MInf=0.7, loc='centers')

### Raccords ####
t = X.connectMatch(t, tol=1.e-6, dim=2)
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           translation=[30.,0.,0.],dim=2)
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           translation=[0.,20,0.],dim=2)



### Ajout des ghost-cells ###
C.addState2Node__(t, 'EquationDimension', 2)
t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
t = T.addkplane(t)
t = T.contract(t, (0,0,0), (1,0,0), (0,1,0),0.025)
t = T.makeDirect(t)

### Construction du tc ###
tc = C.node2Center(t)
tc = X.setInterpData2(t, tc, nature=1, loc='centers', storage='inverse', 
                      sameName=1, method='lagrangian',dim=2)
tc = C.rmVars(tc, 'FlowSolution')
tc = C.rmVars(tc, 'CellN')


# Numerics
numb = {}
numb["temporal_scheme"]    = "explicit"
numb["omp_mode"]           = 0
numb["modulo_verif"]       = 10
numb["rk"]        	   = 3
numb["exp_local"]	   = 2 
numz = {}
numz["time_step"]          = 0.01
numz["scheme"]             = "ausmpred"
numz["niveaux_temps"]      = 1 # Explicit local : nbre niveaux en temps (=1 si explicit global)



internal._setNum2Zones(t, numz) 
internal._setNum2Base(t, numb)

#C.convertPyTree2File(t, 'essai.cgns')

zones = Internal.getNodesFromType2(t, 'Zone_t')

for z in zones:   
      solcenter = Internal.getNodeFromName1(z, 'FlowSolution#Centers') 
      niveau = Internal.getNodeFromName(solcenter, 'niveaux_temps')[1][0][0][0]
      dtloc = Internal.getNodeFromName1(z, '.Solver#define')  # noeud
      level = Internal.getNodeFromName1(dtloc, 'niveaux_temps')  # noeud
      Internal.setValue(level,int(niveau))

(t,tc,metrics) = FastS.warmup(t,tc)

nit=100; times = 0.
timeStep = numz['time_step']
for it in range(nit):
    FastS._compute(t, metrics, it,tc,layer='c')
    if it%10 == 0:
        print('- %d - %f'%(it, times))
        FastS.display_temporal_criteria(t, metrics, it)
    times += timeStep

test.testT(t, 1)
