# Entree: 
# - un maillage 2D (un seul plan) multibloc
# - ou un maillage 3D multibloc 
# avec BCs + raccords + (solution initiale et/ou Reference State)
# Sortie: t.cgns et tc.cgns pour compute.py

import sys
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import numpy
import FastS.PyTree as FastS
import Fast.PyTree  as Fast
import Connector.PyTree as X

Fast.FastC.MX_SYNCHRO= 100000
Fast.FastC.MX_OMP_SIZE_INT= 10000

### Nombre de processus MPI pour la simulation
NP = 3

### Fichier d entree : maillage + CI + CL ###

FILE = 'maillage_lamb.cgns'
#FILE = 'essai.cgns'

t = C.convertFile2PyTree(FILE)

# Get dim
dim = 3
node = Internal.getNodeFromName(t, 'EquationDimension')
if node is not None: dim = Internal.getValue(node)

print('dim= ',dim)

# Solution initiale
eqs = Internal.getNodeFromType(t, 'GoverningEquations_t')
#Model = 'Euler'
if eqs is not None: Model = Internal.getValue(eqs)


### Ajout des ghost-cells pour permettre le calcul de la CFL
C.addState2Node__(t, 'EquationDimension', dim)
t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
if (dim == 2): 
    t = T.addkplane(t)
    t = T.contract(t, (0,0,0), (1,0,0), (0,1,0),0.025)
    t = T.makeDirect(t)



### Met la CFL en chaque cellule dans l arbre t ###
t, exposant_max = FastS.computeCFL_dtlocal(t)
###################################################



### Decoupe le maillage en zones de niveaux en temps different ###
t = FastS._decoupe2(t, exposant_max, NP = NP)
##################################################################


#### On remet les GhostCells car supprimees dans decoupe2
t = Internal.addGhostCells(t, t, 2, adaptBCs=1, fillCorner=0)
if (dim == 2): 
    t = T.addkplane(t)
    t = T.contract(t, (0,0,0), (1,0,0), (0,1,0),0.025)
    t = T.makeDirect(t)


tc = C.node2Center(t)

tc = X.setInterpData3(t, tc, nature=1, loc='centers', storage='inverse', 
                      sameName=1, method='lagrangian',dim=dim)
tc = C.rmVars(tc, 'FlowSolution')
tc = C.rmVars(tc, 'CellN')


C.convertPyTree2File(t, 't.cgns')
C.convertPyTree2File(tc, 'tc.cgns')




