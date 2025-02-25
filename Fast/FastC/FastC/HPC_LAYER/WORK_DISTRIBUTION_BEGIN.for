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

      thread_topology(1) =topo_s(1)
      thread_topology(2) =topo_s(2)
      thread_topology(3) =topo_s(3)

      !cible taille cache bloc applique a la souszone omp
      size_cache(1)       = param_int(CACHEBLCKI)
      size_cache(2)       = param_int(CACHEBLCKJ)
      size_cache(3)       = param_int(CACHEBLCKK)

      do i=1,3
        if(thread_topology(i).ne.1) then
          l            =(ind_dm_omp(2*i)-ind_dm_omp(2*i-1)+1)/2
          size_cache(i)=min( l, size_cache(i) )
          size_cache(i)=max( 1, size_cache(i) )
        endif
      enddo

      thread_pos(3) = 1 + (ithread-1)/(topo_s(1)*topo_s(2))
      l             = ithread -(thread_pos(3)-1)*topo_s(1)*topo_s(2)
      thread_pos(2) = 1 + (l-1)/topo_s(1)
      thread_pos(1) = l - (thread_pos(2)-1)*topo_s(1)
      !on determine synchro_thrread + nbre sous-domaine )cache !bloc)

#if defined(__INTEL_COMPILER)
#if not defined(__INTEL_LLVM_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: crsdm_scater
#endif
#endif
      !in : ind_dm_zone: taille zone, ind_dm_omp: sousdomaine omp
      !out: loop cache bloc (ijkv_sdm), synchro thread
       call crsdm_scater( ndo, topo_s, size_cache,
     &                    synchro_receive_th, synchro_send_th,
     &                    ind_dm_zone   , ind_dm_omp, ijkv_sdm )

c        if(ithread.eq.param_int( IO_THREAD).and.nitrun.eq.0)then
c           write(*,'(a,3i4)')'thread_pos =',thread_pos
c           !write(*,'(a,3i4)')'IJKV socket=',ijkv_thread
c           write(*,'(a,3i4)')'IJKV thread=',ijkv_sdm
c           write(*,'(a,9i4)')'topo thread',thread_topology,ind_dm_omp
c           write(*,'(a,9i4)')'synchro_receive_th',synchro_receive_th
c           write(*,'(a,9i4)')'synchro_send_th',synchro_send_th
c       endif

ccc#include   "../FastC/FastC/HPC_LAYER/topo_cache.for"
ccc#include   "../FastC/FastC/HPC_LAYER/loopBloc_begin.for"
ccc#include   "../FastC/FastC/HPC_LAYER/loop_scater.for"
#include   "../FastC/FastC/HPC_LAYER/verif_loksize.for"

