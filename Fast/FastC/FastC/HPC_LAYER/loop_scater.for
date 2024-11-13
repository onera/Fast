           !on dertemine synchro_thrread + nbre sous-domaine )cache !bloc)

#if defined(__INTEL_COMPILER)
#if not defined(__INTEL_LLVM_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: crsdm_scater
#endif
#endif
           !in : ind_dm_socket: taille zone, ind_dm_omp: sousdomaine omp
           !out: loop cache bloc (ijkv_sdm), synchro thread
           call crsdm_scater( ndo, topo_s, size_cache,
     &                       synchro_receive_th, synchro_send_th,
     &                       ind_dm_zone   , ind_dm_omp, ijkv_sdm )

        if(ithread.eq.param_int( IO_THREAD).and.nitrun.eq.0)then
       !if(ithread.eq.param_int( IO_THREAD).and.mod(nitrun,5).eq.0)then
           !write(*,'(a,4i4)')'ijkboc =',ibloc,jbloc,kbloc,ithread_sock
           !write(*,'(a,4i4)')'ijkboc =',ibloc,jbloc,kbloc
           !write(*,'(a,3i4)')'thread_pos_tmp =',thread_pos_tmp
           !write(*,'(a,8i4)')'lgo    =',lgo,lwait,skip,shift
           !write(*,'(a,9i4)')'topo thread',thread_topology,ind_dm_omp
          if(ibloc*jbloc*kbloc.le.1) then
           write(*,'(a,3i4)')'thread_pos =',thread_pos
           write(*,'(a,3i4)')'IJKV socket=',ijkv_thread
           write(*,'(a,3i4)')'IJKV thread=',ijkv_sdm
           !write(*,'(a,9i4)')'topo socket',socket_topology,ind_dm_socket
           write(*,'(a,9i4)')'topo thread',thread_topology,ind_dm_omp
           write(*,'(a,9i4)')'synchro_receive_th',synchro_receive_th
           write(*,'(a,9i4)')'synchro_send_th',synchro_send_th
          endif
       endif
