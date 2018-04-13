        if(omp_mode.eq.0) then
!DIR$ ATTRIBUTES FORCEINLINE :: indice_boucle_scater
           call indice_boucle_scater(ndo,ibloc,jbloc,kbloc,
     &                            ithread, param_int( IO_THREAD),
     &                            topo_s, size_thread, ijkv_thread,
     &                            ijkvloc, thread_pos_tmp,ind_dm_socket,
     &                            ind_dm_omp )
        else
          ind_dm_omp(1:6)= inddm_omp(1:6)
        endif


           !mise en place du cache bloking + synchro
           if(param_int(KFLUDOM).ne.7)then
           do i=1,3
             if(thread_topology(i)*socket_topology(i).eq.1) then
             !if(thread_topology(i).eq.1) then  
               size_cache(i)=  cache(i)                             
             else
               size_cache(i)=(ind_dm_omp(2*i)-ind_dm_omp(2*i-1)+1)/2
               size_cache(i)=min( cache(i),size_cache(i) )
               size_cache(i)=max(1,size_cache(i))
             endif
           enddo
           else
           size_cache=  cache   
           endif

           !on dertemine synchro_socket

!DIR$ ATTRIBUTES FORCEINLINE :: crsdm
           call crsdm( ndo, size_cache,
     &                 synchro_receive_sock, synchro_send_sock,
     &                 ind_dm_zone, ind_dm_socket, ijkv_sdm ) 

           !on dertemine synchro_thrread + nbre sous-domaine )cache !bloc)

!DIR$ ATTRIBUTES FORCEINLINE :: crsdm_scater
           call crsdm_scater( ndo, ibloc , jbloc , kbloc, topo_s, 
     &                       size_cache,
     &                       synchro_receive_th, synchro_send_th,
     &                       ind_dm_socket     , ind_dm_omp, ijkv_sdm )

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
           write(*,'(a,9i4)')'topo socket',socket_topology,ind_dm_socket
           write(*,'(a,9i4)')'topo thread',thread_topology,ind_dm_omp
           write(*,'(a,9i4)')'synchro_receive_th',synchro_receive_th
           write(*,'(a,9i4)')'synchro_send_th',synchro_send_th
          endif
       endif
