# - compute (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree as Fast
import Fast.Utils as Utils
import FastC.PyTree as FastC
import Initiator.Adim as Adim
import Post.PyTree as P
import Transform.PyTree as T
import numpy
import sys

DIRECTORY_CGNS='cgns_files/'
DIRECTORY_DAT='dat_files/'
DIRECTORY_PLT='plt_files/'

explicit_select = 'explicit'

NP                         = 0     # mpi cores
NIT                        = 1000  # number of iterations
display_probe_freq         = 100   # iteration frequency to display modulo_verif

if explicit_select == 'explicit_local':
    filet    = DIRECTORY_CGNS+'t_lts.cgns'
    filetc   = DIRECTORY_CGNS+'tc_lts.cgns'
else:
    filet    = DIRECTORY_CGNS+'t.cgns'
    filetc   = DIRECTORY_CGNS+'tc.cgns'


restart = False
single_vs_multiple='single' # 'multiple'

NP = Utils.getArgs1()    
if NP > 0: 
    import Converter.Mpi as Cmpi
    import FastS.Mpi as FastS
    rank = Cmpi.rank
    isMPI = True
else: 
    import FastS.PyTree as FastS
    rank = 0
    isMPI = False
if restart:
    if explicit_select == 'explicit_local':
        filet    = DIRECTORY_CGNS+'restart_lts.cgns'
    else:
        filet    = DIRECTORY_CGNS+'restart.cgns'

print('++++ pre-load ++++')
t        = Fast.loadTree(filet , split=single_vs_multiple, mpirun=isMPI)
ts       = None
graph    = None
if isMPI:
    tc,graph = Fast.loadTree(filetc, split='single', graph=True, mpirun=True)
else:
    tc = Fast.loadTree(filetc, split='single')
print('++++ post-load ++++')

it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)

# Numerics
numb = {}
numb["modulo_verif"]       = display_probe_freq
numz = {}
numz["scheme"]             = "ausmpred"

if explicit_select == 'explicit_local':
    numb["temporal_scheme"]    = "explicit_local"
    dt= 5.e-6
    running_cfl                = 0.48
    ## Determines the maximum dt for the given flow field
    dt_max = FastS.set_dt_lts(t,running_cfl,dt)
    numz["time_step"]          = dt_max
else:
    numb["temporal_scheme"]    = "explicit"
    numz["time_step"]          = 5.e-6 # CFL 1

# Numerics
FastC._setNum2Base(t, numb)
## Create for CGNSBase_t
## .Solver#define

FastC._setNum2Zones(t, numz)
## Create for each Zone_t
## .Solver#define

print('++++ pre-warmup ++++')
(t, tc, metrics) = FastS.warmup(t, tc, graph) 
print('++++ post-warmup ++++')

time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)

for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph,layer="Python")
    if it%display_probe_freq == 0:
        if rank == 0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it, format='double')
    time0 += time_step

# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)
if explicit_select == 'explicit_local':
    Fast.saveTree(t ,fileName=DIRECTORY_CGNS+'restart_lts.cgns', split=single_vs_multiple, mpirun=isMPI )#check these files
else:
    Fast.saveTree(t ,fileName=DIRECTORY_CGNS+'restart.cgns', split=single_vs_multiple, mpirun=isMPI )#check these files

