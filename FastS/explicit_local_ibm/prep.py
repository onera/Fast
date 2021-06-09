# IBM preprocessing on Cartesian grids for FastS
import Connector.ToolboxIBM as IBM
import Converter.Internal as Internal
import Converter.PyTree as C
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Fast.PyTree as Fast
import Fast.Utils as Utils
import Generator.PyTree as G
import Initiator.PyTree as I
import Transform.PyTree as T
import sys

sys.path.insert(1, '/stck/ajost/Cassiopee_sator/Apps/PModules/FastS/')
import local_functionalities as LF
DIRECTORY_CGNS='cgns_files/'
DIRECTORY_DAT='dat_files/'
DIRECTORY_PLT='plt_files/'


# Ce script est sequentiel. Si NP > 0, on prepare pour un
# calcul parallele sur NP procs
NP = Utils.getArgs1() # Mpi NP > 0
if NP > 0: print('Preparing for %d processes.'%NP)
if NP > 0: import Distributor2.PyTree as D2
split = 'single' # 'multiple': one file per proc
#=======================================================
# Input data
#=======================================================
# vmin: nb of pts per Cartesian grid
# snears: list of minimum spacing around each surface or
# one value for all surfaces
# dfar: extension of the mesh
# DEPTH: nb of interpolated and ghost cells
vmin = 91; snears = 0.0075; dfar = 5.
DEPTH = 2
# list of body surfaces as separated CGNS bases
bodySurfaceFile = DIRECTORY_CGNS+'naca0012_geom.cgns'
# list of refinement zones as separated CGNS bases
# snearsf size is ensured inside the given surface(s)
refinementSurfaceFile = None # 'refinementBody.cgns'
snearsf = None # 2*snears 
check = False # True: some intermediate files are written for checks
#=======================================================
# End of input data
#=======================================================
#-------------------------------------------------------
# Read body mesh
tb = C.convertFile2PyTree(bodySurfaceFile)

#Set new equation dimension, reference state, governing equations, flow solution, flowequation set, etc.
Mus = 1.78938e-5*3.866
LF.set_new_control_vars(tb,uinf=60.0,angle=0,isTurb=False,isEuler=True,dimPb = 2,Mus=Mus)

#-------------------------------------------------------
# Refinement surfaces in the fluid
#-------------------------------------------------------
try: tbox = C.convertFile2PyTree(refinementSurfaceFile)
except: tbox=None # no refinement surface

#--------------------------------------------------------
# Get Reference State and model from body pyTree
model = Internal.getNodeFromName(tb, 'GoverningEquations')
if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
# model: Euler, NSLaminar, NSTurbulent
model = Internal.getValue(model)

# reference state
refstate = C.getState(tb)
# dimension du pb
dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
dimPb = Internal.getValue(dimPb)
if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

#--------------------------------------------------------
# Generates the full Cartesian mesh
t = IBM.generateIBMMesh(tb,vmin,snears,dfar,DEPTH=DEPTH,tbox=tbox,snearsf=snearsf,check=check)

#------------------------------------------------
# Add reference state to the pyTree and init fields
# Add viscosity if model is not Euler
C._addState(t, state=refstate)
C._addState(t, 'GoverningEquations', model)
C._addState(t, 'EquationDimension', dimPb)
if check: C.convertPyTree2File(t, DIRECTORY_CGNS+'mesh1.cgns')

#----------------------------------------
# Computes distance field
#----------------------------------------
if dimPb == 2:
    z0 = Internal.getNodeFromType2(t,"Zone_t")
    bb = G.bbox(z0); dz = bb[5]-bb[2]
    tb2 = C.initVars(tb,'CoordinateZ',dz*0.5)
    DTW._distance2Walls(t,tb2,type='ortho',signed=0, dim=dimPb,loc='centers')
else:
    DTW._distance2Walls(t,tb,type='ortho',signed=0, dim=dimPb,loc='centers')
    
#----------------------------------------
# Create IBM info
#----------------------------------------
t,tc = IBM.prepareIBMData(t, tb, frontType=1)


#------------------------------------------------------
# distribute the mesh over NP processors
if NP > 0:
    print('distribution over %d processors'%NP)
    stats = D2._distribute(t, NP)
    if check: print(stats)

#----------------------------------------
# arbre donneur
#----------------------------------------
D2._copyDistribution(tc,t)
Fast.save(tc, DIRECTORY_CGNS+'tc.cgns', split='single', NP=-NP)

#----------------------------------------
# Extraction des coordonnees des pts IBM
#----------------------------------------
if check:
    tibm = IBM.extractIBMInfo(tc)
    C.convertPyTree2File(tibm, DIRECTORY_CGNS+'IBMInfo.cgns')
    tibm = None

#----------------------------------------
# arbre de calcul 
#----------------------------------------
del tc # to free memory before initialization
I._initConst(t, loc='centers')
if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)

Fast.save(t, DIRECTORY_CGNS+'t.cgns', split=split, NP=-NP)
