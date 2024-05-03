// Les protos fortran
extern "C"
{
  void nature_geom_dom_( E_Int* nijk_xyz  ,
                         E_Int& ndimdx_xyz, E_Float* x         , E_Float* y       , E_Float* z       , E_Int* degener   ,
                         E_Int& lale      , E_Float* param_real, E_Int* nijk_mtr  , E_Int& ndimdx_mtr,
                         E_Int& neq_ij    , E_Int& neq_k     , E_Int& typ_zone);

  void skmtr_( E_Int& nd      , E_Int* ipt_param_int, E_Float*  ipt_param_real,
               E_Float* x     , E_Float* y          , E_Float* z  , E_Int* degener   , E_Float* dist,
               E_Float* ti    , E_Float* tj         , E_Float* tk , E_Float* ti0     , E_Float*  tj0, E_Float* tk0, E_Float* vol, E_Float* venti, E_Float* ventj, E_Float* ventk,
               E_Int* ijkv_sdm  , 
               E_Int* ind_sdm   , E_Int* ind_coe          , E_Int* ind_grad    , E_Int* ind_dm_zone,  E_Int* ind_dm_socket,  E_Int* ind_dm_omp,  E_Int* socket_topology,
               E_Int& ithread_sock, E_Int& thread_parsock, E_Int& Nbre_thread_actif, E_Int&  Nbre_socket, E_Int& socket     ,
               E_Int& ithread);

   void     cp_tijk_( E_Int* param_int, E_Float* iptx , E_Float* ipty , E_Float* iptz,
                      E_Float* ipti   , E_Float* iptj , E_Float* iptk ,
                      E_Float* iptx0  , E_Float* ipty0, E_Float* iptz0, E_Int* ind_loop );

   void     cp_vol_( E_Int* param_int, E_Float* iptx , E_Float* ipty , E_Float* iptz,
                     E_Float* ipti   , E_Float* iptj , E_Float* iptk ,
                     E_Float* iptx0  , E_Float* ipty0, E_Float* iptz0, E_Float* iptvol,  E_Int* ind_loop );

   void     tijk_extrap_( E_Int& ndimdx_mtr, E_Int& ndimdx_xyz, E_Int* nijk_xyz,  E_Int* nijk_mtr, E_Int& neq_ij , E_Int& neq_k,
                          E_Int* ind_zone , E_Int*  degen   , 
                          E_Float* ipti   , E_Float* iptj , E_Float* iptk ,
                          E_Float* iptx0  , E_Float* ipty0, E_Float* iptz0,  E_Float* iptvol);

   void     dist_extrap_( E_Int& ndimdx   , E_Int& ndimdx_xyz, E_Int* nijk,  E_Int* nijk_xyz,
                          E_Int* ind_zone , E_Int*  degen    , E_Float* iptdist);

  void indice_boucle_lu_(E_Int& ndom , E_Int& ndsdm           , E_Int& nsdom_lu, E_Int& icp ,
                        E_Int* ind_dm, E_Int*  thread_topology, E_Int*  ind_sdm);

  void copynuma_( E_Int*  ind_loop,  E_Int& ni,  E_Int& nj, E_Int& shift, E_Int& ific, E_Int& jfic, E_Int& kfic,
                 E_Float* target       , E_Float* src );

  void cpmys_rij_( E_Int& ndo           , E_Int& ndimdx_my , E_Int& neq_my    , E_Int& neq_grad   ,
                   E_Int& lthermique    , E_Int& ltensrey  , E_Int& lcyl,
                   E_Int* ipt_nijk_my   , E_Int* iptmoy_param  ,  E_Int* param_int,
                   E_Int* ipt_ind_dm_omp,
                   E_Float&  gamma      , E_Float&  cv         ,  E_Float&  prandtl  ,
                   E_Float* iptro       ,
                   E_Float* iptmut      , 
                   E_Float* iptx        , E_Float* ipty        , E_Float* iptz           ,
                   E_Float* ipti        , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                   E_Float* ipti_df     , E_Float* iptj_df     , E_Float* iptk_df        , E_Float* iptvol_df      , 
                   E_Float* iptromoy    );

  void verrou_c_(E_Int* lok , E_Int& type );

  void flush_real_(  E_Int& size, E_Float* tab);

  void indice_cl_sdm_( E_Int& dir    , E_Int& npass  , E_Int& lskip , E_Int& typbc,
                       E_Int& ific   , E_Int& kfic   , 
                       E_Int* ijkv   , E_Int* ind_fen, E_Int* ind_dm, E_Int* ind_dm_thread, 
                       E_Int* ind_CL , E_Int* ind_CL119);

  void init_ssiter_bloc_(E_Int& nd               , E_Int& nitcfg     , E_Int& nssiter , E_Int& nitrun ,
                         E_Int&  lssiter_loc     , E_Int& itypcp     , E_Int& flag_res,
                         E_Int* ijkv             , E_Int*  ijkv_lu   , E_Int*  ijk_lu , E_Int*  size_ssdom,
                         E_Int& mx_ssdom_lu      , E_Int* iskip_lu   , E_Int* iptdtloc,
                         E_Int*   ipt_ind_dm     , E_Int*   ipt_nidom_loc      , E_Int& it_bloc             , E_Int*   ipt_nisdom_residu,
                         E_Int*   ipt_it_lu_ssdom, E_Int*   ipt_it_target_ssdom, E_Int*   ipt_it_target_old , E_Int*   ipt_no_lu, E_Int*   param_int );

 }
