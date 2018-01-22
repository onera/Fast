// Les protos fortran
extern "C"
{
  void nature_geom_dom_( E_Int* nijk_xyz  ,
                         E_Int& ndimdx_xyz, E_Float* x         , E_Float* y       , E_Float* z       , E_Int* degener   ,
                         E_Int& lale      , E_Float* param_real, E_Int* nijk_mtr  , E_Int& ndimdx_mtr,
                         E_Int& neq_ij    , E_Int& neq_k     , E_Int& typ_zone);

  void skmtr_( E_Int& nd      , E_Int* ipt_param_int, E_Float*  ipt_param_real,  E_Float* rot_ale,
               E_Float* x     , E_Float* y          , E_Float* z  , E_Int* degener   ,
               E_Float* ti    , E_Float* tj         , E_Float* tk , E_Float* ti0     , E_Float* tj0, E_Float* tk0, E_Float* vol, E_Float* venti, E_Float* ventj, E_Float* ventk,
               E_Int* ijkv_sdm  , 
               E_Int* ind_sdm   , E_Int* ind_coe          , E_Int* ind_grad    , E_Int* ind_dm_zone,  E_Int* ind_dm_socket,  E_Int* ind_dm_omp,  E_Int* socket_topology,
               E_Int& ithread_sock, E_Int& thread_parsock, E_Int& Nbre_thread_actif, E_Int&  Nbre_socket, E_Int& socket     ,
               E_Int& ithread);

  void skmtr_df_( E_Int* nijk_xyz   , E_Int* nijk_mtr , 
                  E_Int& ndimdx_xyz , E_Float* x      , E_Float* y     , E_Float* z      ,
                  E_Float* tmpi     , E_Float* tmpj   , E_Float* tmpk  , E_Float* tmpi2  , E_Float* tmpj2, E_Float* tmpk2 , E_Float* mtr,
                  E_Int& neq_ij     , E_Int& neq_k    ,
                  E_Int& ndimdx_mtr , E_Float* ti     , E_Float* tj    , E_Float* tk     , E_Float* vol);


 // void init_ssiter_bloc_(E_Int& nd    , E_Int&  lssiter_loc      );

  void init_ssiter_bloc_(E_Int& nd               , E_Int& nitcfg     , E_Int& nssiter , 
                         E_Int&  lssiter_loc     , E_Int& itypcp     ,
                         E_Int* ijkv             , E_Int*  ijkv_lu   , E_Int*  ijk_lu , E_Int*  size_ssdom,
                         E_Int& mx_ssdom_lu      , E_Int* iskip_lu   , 
                         E_Int*   ipt_ind_dm     , E_Int*   ipt_nidom_loc      , E_Int& it_bloc             , E_Int*   ipt_nisdom_residu,
                         E_Int*   ipt_it_lu_ssdom, E_Int*   ipt_it_target_ssdom, E_Int*   ipt_it_target_old , E_Int*   ipt_no_lu, E_Int*   param_int );

  void cprdu3s1_(E_Int& nd               , E_Int& nitcfg     , E_Int& nssiter            , E_Int& neq              ,  E_Int& ndimdx       , E_Int& ndim_rdm  , E_Int& nitrun, 
                 E_Int& iflw             , E_Int& iles       , E_Int& ithread            , E_Int& omp_mode         ,  E_Int& mx_ssdom_lu  ,
                 E_Int& it_bloc          , E_Float& epsi     , E_Int*  size_ssdom        , E_Int* ipt_nisdom_residu,
                 E_Int*  nijk            , E_Int* ijkv       , E_Int*  ijkv_lu           , 
                 E_Int* ind_loop         ,
                 E_Int* ipt_it_lu_ssdom  , E_Int* ipt_it_target_ssdom,  E_Int* ipt_it_target_old, E_Int* ipt_it_temp_ssdom,  E_Int* ipt_no_lu, 
                 E_Float* iptdrodm       , E_Float* rdm  );


  void indice_boucle_lu_(E_Int& ndom , E_Int& ndsdm           , E_Int& nsdom_lu, E_Int& icp ,
                        E_Int* ind_dm, E_Int*  thread_topology, E_Int*  ind_sdm);

  void copy_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,   E_Float* iptro, E_Float* stock, E_Int& ind, E_Int& nzone);

  
  void interpolation_(E_Int& idir , E_Int* ipt_param_int , E_Float* ipt_param_real ,E_Int* ind_loop, E_Float* iptro_tmp , E_Float* iptro);

   
  void copynuma_( E_Int*  ind_loop,  E_Int& ni,  E_Int& nj, E_Int& shift, E_Int& ific, E_Int& jfic, E_Int& kfic,
                 E_Float* target       , E_Float* src );

  void copyflux_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,   E_Float*drodm, E_Float* stock, E_Int& ind, E_Int& nzone);

  void init_ventijk_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   ,
                      E_Int* ipt_param_int      , E_Float* ipt_param_real  , 
                      E_Int* ipt_ijkv_sdm       ,
                      E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
                      E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                      E_Float* ipti_df          , E_Float* iptj_df     , E_Float* iptk_df        , 
                      E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                      E_Float* iptx             , E_Float* ipty         , E_Float* iptz    );


  void post_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& neq_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real  ,  E_Float& tke      , E_Float& enst     ,  E_Int& compteur  ,
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptgrad   );

  void post_grad_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& neq_grad, E_Int& order,
              E_Int* ipt_param_int      , E_Float* ipt_param_real  ,
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptgrad   );


  void post_q_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ,  E_Float* iptenst ,  E_Float* iptvort   );

  void post_qprime_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            , E_Float* iptro_m1    ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ   );

  void post_qprime_1rot_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif , E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   ,
              E_Int& order              ,  E_Int& dim_grad          , E_Int& var1, E_Int& var2       ,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            , E_Float* iptro_m1    ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ , E_Float* iptrot );


  void post_enst_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ,  E_Float* iptenst ,  E_Float* iptvort   );

  void post_q_enst_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ,  E_Float* iptenst ,  E_Float* iptvort   );


  void post_drodt_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            , E_Float*rop_n         , E_Float* rop_n1         , E_Float*drodt );


  void viles_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& neq_grad, 
              E_Int* ipt_param_int      , E_Float* ipt_param_real  ,
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol          ,  E_Float* iptmut, E_Float* iptdist, E_Float* iptrot   );

  void navier_stokes_struct_( E_Int& ndo    , E_Int& nidom            , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , 
                              E_Int& lssiter_verif  , E_Int& nptpsi            , E_Int& nitcfg  , E_Int& nitrun     , E_Int& first_it   ,  E_Int& nb_pulse  , E_Int&   flagCellN  ,
                              E_Int* ipt_param_int  , E_Float* ipt_param_real  ,
                              E_Float& temps        , E_Int* ipt_tot,   
                              E_Int* ipt_ijkv_sdm       ,
                              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
                              E_Float* ipt_cfl        ,
                              E_Float* iptx             , E_Float* ipty         , E_Float* iptz          , E_Float* iptCellN       ,
                              E_Float* iptro            , E_Float* iptro_m1    , E_Float* iptrotmp       , E_Float* iptro_ssiter   ,
                              E_Float* iptmut           , E_Float* iptdist     ,
                              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                              E_Float* ipti_df          , E_Float* iptj_df     , E_Float* iptk_df        , E_Float* iptvol_df      , 
                              E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                              E_Float* iptwig           , E_Float* iptstat_wig , E_Float* iptrot         ,
                              E_Float* iptdrodm         , E_Float* iptcoe      , E_Float* iptdelta       , E_Float* iptfd          , E_Float* iptro_zgris    , E_Float* iptro_res );

  void invlu_(                E_Int& ndo      , E_Int& nitcfg      , E_Int& nitrun   , E_Int*  param_int , E_Float* param_real,
                              E_Int* ipt_sdm            ,
                              E_Float* iptrotmp         , E_Float* iptro_ssiter   ,
                              E_Float* iptdrodm         , 
                              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , 
                              E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                              E_Float* iptcoe    );




   void     bvbs_extrapolate_( E_Int& idir , E_Int& lrhs    , E_Int& eq_deb     ,E_Int* param_int , E_Int* ind_loop  , E_Float& nutildeinf,  E_Float* iptro);

   void     bvbs_periodique_( E_Int& idir  , E_Int& lrhs    , E_Int* param_int  , E_Int* ind_loop, E_Float* iptro);

   void     bvbs_periodique_azimuthal_( E_Int& idir  , E_Int& lrhs    , E_Int* param_int  , E_Int* ind_loop, E_Float* iptro, E_Float* iptdata);

   void     bvbs_wall_inviscid_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                 E_Float* iptventi  , E_Float* iptijk  , E_Float* iptro);

   void     bvbs_wall_viscous_adia_( E_Int& idir      , E_Int& lrhs      , E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                     E_Float* iptventi, E_Float* iptijk  , E_Float* iptro);

   void     bvbs_wall_viscous_transition_( E_Int& idir      , E_Int& lrhs      ,   E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                           E_Float* param_real,
                                           E_Float* iptx    , E_Float* ipty     , E_Float* iptz  ,
                                           E_Float* iptventi, E_Float* iptijk   , E_Float* iptro);


   void     bvbs_inflow_supersonic_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                     E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                     E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);

   void     bvbs_farfield_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                            E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                            E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);

   void     bvbs_outflow_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                            E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                            E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);

   void     bvbs_inflow_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                           E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                           E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);


   void     bvbs_outpres_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr,  E_Int* param_int ,E_Int* ind_loop  ,
                            E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                            E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* data_pres,
			    E_Int& size_data, E_Int& inc_bc);

//   void     bvbs_updatepressure_(  E_Int& idir  , E_Int& ithread, E_Int* ind_avg, E_Int* ind_mjr, E_Int* param_in ,E_Int* ind_loop ,
//                            E_Float* param_real, E_Float* ipty   , E_Float* iptz  , E_Float* iptro, E_Float* data_pres,
//			    E_Int& size_data, E_Float* vteta, E_Float* roteta, E_Int& inc_bc);

   void     bvbs_inflow_newton_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                  E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                  E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, 
                                  E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	          E_Int& size_data   , E_Int& inc_bc     , E_Int& size_work);

   void     bvbs_inflow_fich_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, 
                                E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	        E_Int& size_data   , E_Int& inc_bc     , E_Int& size_work);


   void indice_cl_sdm_( E_Int& dir    , E_Int& npass  , E_Int& lskip , E_Int& typbc,
                        E_Int& ific   , E_Int& kfic   , 
                        E_Int* ijkv   , E_Int* ind_fen, E_Int* ind_dm, E_Int* ind_dm_thread, 
                        E_Int* ind_CL , E_Int* ind_CL119);

   void     correct_coins_( E_Int& ndo ,  E_Int* param_int, E_Int* ind_loop, E_Float* iptro);


   void     move_domx_( E_Int& ndo ,  E_Int* param_int, E_Float* param_real, E_Int* ind_loop, 
                        E_Float* iptx  , E_Float* ipty , E_Float* iptz, 
                        E_Float* iptx0 , E_Float* ipty0, E_Float* iptz0  );

   void     mjr_ale_3dhomocart_( E_Int& ndo  , E_Int* param_int   , E_Float* param_real,
                        E_Int& socket        , E_Int& Nbre_socket , E_Int& thread_sock , 
                        E_Int& thread_parsock, E_Int* ind_socket  , E_Int* topo,
                        E_Float* iptx  , E_Float* ipty , E_Float* iptz, 
                        E_Float* ipti  , E_Float* iptj , E_Float* iptk ,
                        E_Float* ipti0 , E_Float* iptj0, E_Float* iptk0, E_Float* iptvol ,E_Float* venti, E_Float* ventj, E_Float* ventk);


   void     dpssiter_(  E_Int& nitrun , E_Int& neq , E_Int& nssiter, E_Int& iflw, E_Int& iles, E_Int& lft, char*, E_Int& size_name, E_Float* rdm, E_Float* cvg_ptr);

   void     conv2pytree_(E_Float* cvg_pt, E_Int& nitrun, E_Int& neq, E_Int* LastRec, char* name, E_Int& size_name, E_Int& lft, E_Int& nrec, E_Int& nd, E_Int* Itnum, E_Float* Res_L2, E_Float* Res_oo);

  void cpmys_rij_( E_Int& ndo           , E_Int& ndimdx    , E_Int& ndimdx_mtr , E_Int& ndimdx_my  , 
                   E_Int& neq           , E_Int& neq_my    , E_Int& neq_grad   , E_Int& neq_ij     , E_Int& neq_k  , E_Int& imtr       ,
                   E_Int& lthermique    , E_Int& ltensrey  , 
                   E_Int* ipt_nijk_mtr  , E_Int* ipt_nijk  , E_Int* ipt_nijk_my, E_Int* ipt_ijkv   ,
                   E_Int* iptmoy_param  ,
                   E_Int* ipt_ind_dm_omp,
                   E_Float&  gamma      , E_Float&  cv         ,  E_Float&  prandtl  ,
                   E_Float* iptro       ,
                   E_Float* iptmut      , 
                   E_Float* ipti        , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                   E_Float* ipti_df     , E_Float* iptj_df     , E_Float* iptk_df        , E_Float* iptvol_df      , 
                   E_Float* iptromoy    );

  void bceffort_( E_Int& ndo        ,  E_Int& ithread     ,  E_Int* param_int, E_Float* param_real,  E_Int* param_int_eff,
                  E_Int* ind_loop   ,  E_Float* effort_omp,  E_Float* xyz_ref,
                  E_Float*  iptro   ,  E_Float* iptflu    ,  E_Float* iptwig , E_Float* iptmut,   
                  E_Float* iptx     ,  E_Float* ipty      ,  E_Float* iptz   ,
                  E_Float* ipti     ,  E_Float* iptj      ,  E_Float* iptk   , E_Float* iptvol,
                  E_Float* venti    ,  E_Float* ventj     ,  E_Float* ventk );

  void sfd_( E_Int* param_int, E_Float* param_real, E_Int& nitrun, E_Int* ind_loop, E_Float* ipt_CL, E_Float* iptrof, E_Float* iptcoe, E_Float* iptvol );

  void dpj_dpw_calc_meth_( E_Int& ndo        ,  E_Int& ithread       ,  E_Int* param_int, E_Float* param_real,  E_Int* param_int_eff,
                  E_Int* ind_loop      ,  E_Float* xyz_ref     ,
                  E_Float*  iptro      ,  E_Float* iptwig      , E_Float* iptmut  ,   
                  E_Float* iptx        ,  E_Float* ipty        ,  E_Float* iptz   ,
                  E_Float* dpCdp_dpW   ,  E_Float* dpClp_dpW   ,
                  E_Float* ipti        ,  E_Float* iptj        ,  E_Float* iptk   , E_Float* iptvol ,
                  E_Float* venti       ,  E_Float* ventj       ,  E_Float* ventk  , E_Float& cosAoA , E_Float& sinAoA, E_Float& surfinv); 

  void rhs_adjoint_( E_Int& ndo    , E_Int& nidom            , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , 
                              E_Int& lssiter_verif  , E_Int& nptpsi            , E_Int& nitcfg  , E_Int& nitrun     , E_Int& first_it   ,  E_Int& nb_pulse  , E_Int&   flagCellN  ,
                              E_Int* ipt_param_int  , E_Float* ipt_param_real  ,
                              E_Float& temps        , E_Int* ipt_tot,   
                              E_Int* ipt_ijkv_sdm       ,
                              E_Int* ipt_sdm            , E_Int* ipt_coe        , E_Int* ipt_grad         , 
                              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
                              E_Float* ipt_cfl        ,
                              E_Float* iptx             , E_Float* ipty         , E_Float* iptz          , E_Float* iptCellN       ,
                              E_Float* iptro            ,
                              E_Float* iptmut           , E_Float* iptdist     ,
                              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                              E_Float* ipti_df          , E_Float* iptj_df     , E_Float* iptk_df        , E_Float* iptvol_df      , 
                              E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                              E_Float* iptwig           , E_Float* iptstat_wig ,
                              E_Float* iptdrodm         , E_Float* iptcoe     , 
                              E_Float* dpJ_dpW          , E_Float* rhsIter    , E_Float* Adj);
}
