import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree  as Fast
import Fast.Utils as Utils
import FastC.PyTree  as FastC
import FastS.PyTree as FastS
import Transform.PyTree as T
import math
import numpy
import sys

## NOTE::
# Input: 
# - t.cgns & tc.cgns file that conforms to the output of FastS prep.py for IBM
## The latter is needed for the IBM interpolations
# Output:
# - t.cgns & tc.cgns compute.py

## The current file just adds the needed info in the t & tc PyTrees.
## As such a correct t must be provided w/ an initial flow field

DIRECTORY_CGNS='cgns_files/'
DIRECTORY_DAT='dat_files/'
DIRECTORY_PLT='plt_files/'

Fast.FastC.MX_SYNCHRO      = 100000
Fast.FastC.MX_OMP_SIZE_INT = 10000

NP = Utils.getArgs1() 

filet  = DIRECTORY_CGNS+'restart.cgns'
t      = C.convertFile2PyTree(filet)
filetc = DIRECTORY_CGNS+'tc.cgns'
tc     = C.convertFile2PyTree(filetc)

### Met la CFL en chaque cellule dans l arbre t ###
t, exposant_max = FastS.computeCFL_dtlocal(t)
t,tc = FastS._decoupe4(t,tc,exposantMax=exposant_max,isOctree=True)
                    
C.convertPyTree2File(t , DIRECTORY_CGNS+'t_lts.cgns')
C.convertPyTree2File(tc, DIRECTORY_CGNS+'tc_lts.cgns')

 


