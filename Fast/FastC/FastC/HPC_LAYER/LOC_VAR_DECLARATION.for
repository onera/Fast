      logical lerr

      character*7 omp_init,omp_wait,omp_go,omp_wait_lu

      INTEGER_E i, icache,jcache,kcache,lmin,
     & size_max, size_loc,thread_parsock,
     & thread_parsock_actif,extended_range,
     & lok_shap_sock(4),thread_topology(3),
     & lok_shap(4), size_cache(3),
     & socket_pos(3), synchro_receive_sock(3),synchro_send_sock(3),
     & synchro_receive_th(3),synchro_send_th(3),ipt_lok_sock,
     & ithread_sock,ipt_lok,size_max_sock,neq_lok,taille,
     & ijkv_thread(3),kGbloc,jGbloc,iGbloc,ip,jp,kp,lth,
     & ibloc,jbloc,kbloc,ijkvloc(3),skip(3),shift(3),test(3),lwait,lgo,
     & size_thread(3),thread_pos(3),thread_pos_tmp(3),sens(3),
     & size_target(3), cache(3)

      INTEGER_E ind_coe(6),ind_grad(6),ind_sdm(6),ind_rhs(6),ind_mjr(6),
     & ind_ssa(6), ind_hrr(6), ind_gcb(6), ind_src(6)

#if CHECK_SPLIT > 1
      INTEGER_E   inc1,inc2,inc11,inc22
#endif

