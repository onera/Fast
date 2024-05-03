      !cible taille souszone omp si socket!=1
      size_target(1) =  2048
      size_target(2) =  3
      size_target(3) =  7
      
      !cible taille cache bloc applique a la souszone omp
      cache(1)       = param_int(CACHEBLCKI)
      cache(2)       = param_int(CACHEBLCKJ)
      cache(3)       = param_int(CACHEBLCKK)

#if defined(__INTEL_COMPILER)
#if not defined(__INTEL_LLVM_COMPILER)
!DIR$ ATTRIBUTES FORCEINLINE :: topo_scater
#endif
#endif
      call topo_scater(ndo, ithread, socket, lmin,
     &                  param_int(KFLUDOM),
     &                  thread_parsock, thread_parsock_actif,
     &                  ithread_sock, socket_topology,
     &                  size_target, ind_dm_zone,
     &                  topo_s, ijkvloc, thread_pos, socket_pos,
     &                  thread_topology, size_thread, ijkv_thread)

      sens(:)=1
      skip(:)=1

c      if(topo_s(1).ne.1.and.ijkv_thread(1).ne.1
c     &                 .and.mod(ijkv_thread(1),2).eq.2 ) skip(1)=2
c      if(topo_s(2).ne.1.and.ijkv_thread(2).ne.1
c     &                 .and.mod(ijkv_thread(2),2).eq.2 ) skip(2)=2
c      if(topo_s(3).ne.1.and.ijkv_thread(3).ne.1
c     &                 .and.mod(ijkv_thread(3),2).eq.2 ) skip(3)=2

      lok_shap_sock(1)  = 2* min( socket_topology(1)-1, 1 )
      lok_shap_sock(2)  = 2* min( socket_topology(2)-1, 1 )
      lok_shap_sock(3)  = 2* min( socket_topology(3)-1, 1 )

      size_max_sock=ijkv_thread(2)*ijkv_thread(3)*lok_shap_sock(1)
      size_loc     =ijkv_thread(1)*ijkv_thread(3)*lok_shap_sock(2)
      size_max_sock= max(size_max_sock,size_loc)
      size_loc     =ijkv_thread(1)*ijkv_thread(2)*lok_shap_sock(3)

      lok_shap_sock(4)= max(size_max_sock,size_loc)/2

      lok_shap_sock(3)= max(1,lok_shap_sock(1)+lok_shap_sock(2)
     &                                        +lok_shap_sock(3) )

      ipt_lok_sock = Nbre_thread_actif + 1
      ipt_lok      = ipt_lok_sock  + lok_shap_sock(4)*lok_shap_sock(3)
     &                                               *Nbre_socket
