import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree  as Fast
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Transform.PyTree as T
import numpy
import math
import sys

## NOTE::
# Input: 
# - t.cgns file that conforms to the output of FastS prep.py 
# Output:
# - t.cgns & tc.cgns compute.py

## The current file just adds the needed info in the t & tc PyTrees.
## As such a correct t must be provided w/ an initial flow field


DIRECTORY_CGNS='cgns_files/'
DIRECTORY_DAT='dat_files/'
DIRECTORY_PLT='plt_files/'

Fast.FastC.MX_SYNCHRO     = 100000
Fast.FastC.MX_OMP_SIZE_INT= 10000

### Nombre de processus MPI pour la simulation
NP = 0

##Input file
filename = 'restart.cgns'

FILE = DIRECTORY_CGNS+filename
t = C.convertFile2PyTree(FILE)

### Met la CFL en chaque cellule dans l arbre t ###
t, exposant_max = FastS.computeCFL_dtlocal(t)

### Decoupe le maillage en zones de niveaux en temps different ###
t,tc = FastS._decoupe4(t, exposant_max, NP = NP,taille_bloc=25)

C.convertPyTree2File(t , DIRECTORY_CGNS+'t2DLTS.cgns')
C.convertPyTree2File(tc, DIRECTORY_CGNS+'tc2DLTS.cgns')





 


