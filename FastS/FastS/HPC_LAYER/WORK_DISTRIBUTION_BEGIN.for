#define CHECK_SPLIT 0
#define CHECK_BLOCK 1
#if CHECK_SPLIT > 0
      tot(:,ithread)=0
#endif

      extended_range= 0
      omp_init      = 'init   '
      omp_wait      = 'wait   '
      omp_go        = 'go     '
      omp_wait_lu   = 'wait_lu'
      lerr          =.false.

      thread_parsock =  Nbre_thread_actif/Nbre_socket

#include   "FastS/HPC_LAYER/topo_cache.for"

#include   "FastS/HPC_LAYER/loopBloc_begin.for"

#include   "FastS/HPC_LAYER/loop_scater.for"

#include   "FastS/HPC_LAYER/verif_loksize.for"
