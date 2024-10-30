      !cible taille souszone omp si socket!=1
      !size_target(1) =  2048
      !size_target(2) =  3
      !size_target(3) =  7
      
      !cible taille cache bloc applique a la souszone omp
      size_cache(1)       = param_int(CACHEBLCKI)
      size_cache(2)       = param_int(CACHEBLCKJ)
      size_cache(3)       = param_int(CACHEBLCKK)

      thread_pos(3) = 1 + (ithread-1)/(topo_s(1)*topo_s(2))
      l             = ithread -(thread_pos(3)-1)*topo_s(1)*topo_s(2)
      thread_pos(2) = 1 + (l-1)/topo_s(1)
      thread_pos(1) = l - (thread_pos(2)-1)*topo_s(1)


#if defined(__INTEL_COMPILER)
#if not defined(__INTEL_LLVM_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: topo_scater
#endif
#endif
c      !!out: thread_pos, socket_pos, thread_topology, size_thread, ijkv_thread
c      call topo_scater(ndo, ithread, socket, lmin,
c     &                  param_int(KFLUDOM),
c     &                  thread_parsock, thread_parsock_actif,
c     &                  ithread_sock, socket_topology,
c     &                  size_target, ind_dm_zone,
c     &                  topo_s, ijkvloc, thread_pos, socket_pos,
c     &                  thread_topology, size_thread, ijkv_thread)

c      sens(:)=1
c      skip(:)=1

c      lok_shap_sock(1)  = 2* min( socket_topology(1)-1, 1 )
c      lok_shap_sock(2)  = 2* min( socket_topology(2)-1, 1 )
c      lok_shap_sock(3)  = 2* min( socket_topology(3)-1, 1 )

c      size_max_sock=ijkv_thread(2)*ijkv_thread(3)*lok_shap_sock(1)
c      size_loc     =ijkv_thread(1)*ijkv_thread(3)*lok_shap_sock(2)
c      size_max_sock= max(size_max_sock,size_loc)
c      size_loc     =ijkv_thread(1)*ijkv_thread(2)*lok_shap_sock(3)

c      lok_shap_sock(4)= max(size_max_sock,size_loc)/2

c      lok_shap_sock(3)= max(1,lok_shap_sock(1)+lok_shap_sock(2)
c     &                                        +lok_shap_sock(3) )

c      ipt_lok_sock = Nbre_thread_actif + 1
c      ipt_lok      = ipt_lok_sock  + lok_shap_sock(4)*lok_shap_sock(3)
c     &                                               *Nbre_socket
