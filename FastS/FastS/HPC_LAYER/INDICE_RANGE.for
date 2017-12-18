!DIR$ ATTRIBUTES FORCEINLINE :: indice_boucle_ssdom
           call indice_boucle_ssdom(ndo, extended_range,
     &                              ibloc , jbloc , kbloc,
     &                              icache, jcache, kcache,
     &                              topo_s, ithread_sock,thread_pos_tmp,
     &                              size_cache,
     &                              synchro_receive_sock,
     &                              synchro_receive_th  ,
     &                              synchro_send_sock,
     &                              synchro_send_th  ,
     &                              param_int(NIJK), param_int(IJKV),
     &                              ind_dm_zone, ind_dm_socket,
     &                              ind_dm_omp, ijkv_sdm,
     &                              ind_sdm , ind_coe,
     &                              ind_grad, ind_rhs, ind_mjr)

#if CHECK_SPLIT > 0
#include "FastS/HPC_LAYER/check_split0.for"
#endif
