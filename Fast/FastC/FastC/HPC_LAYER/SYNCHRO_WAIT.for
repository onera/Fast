              call synchro_omp_scater(param_int, ithread,
     &                          Nbre_thread_actif,
     &                          lok_shap, neq_lok,
     &                          ithread, thread_topology, thread_pos,
     &                          synchro_receive_th, synchro_send_th,
     &                          icache, jcache, kcache, ijkv_sdm,
     &                          size_cache,
     &                          ind_dm_omp,
     &                          lok(ipt_lok), omp_wait )
