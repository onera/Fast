// Les protos fortran
extern "C"
{
  void k6stretch_(E_Float& t1, E_Float& t2, E_Float* sn, E_Int& nbp,E_Float& dsm, E_Float& dsp, E_Int& ityp,E_Int& inewt);

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

  void skmtr_df_( E_Int* nijk_xyz   , E_Int* nijk_mtr , 
                  E_Int& ndimdx_xyz , E_Float* x      , E_Float* y     , E_Float* z      ,
                  E_Float* tmpi     , E_Float* tmpj   , E_Float* tmpk  , E_Float* tmpi2  , E_Float* tmpj2, E_Float* tmpk2 , E_Float* mtr,
                  E_Int& neq_ij     , E_Int& neq_k    ,
                  E_Int& ndimdx_mtr , E_Float* ti     , E_Float* tj    , E_Float* tk     , E_Float* vol);


  void init_ssiter_bloc2_(E_Int& nd               , E_Int& nitcfg     , E_Int& nssiter , 
                         E_Int&  lssiter_loc     , E_Int& itypcp     ,
                         E_Int* ijkv             , E_Int*  ijkv_lu   , E_Int*  ijk_lu , E_Int*  size_ssdom,
                         E_Int& mx_ssdom_lu      , E_Int* iskip_lu   , 
                         E_Int*   ipt_ind_dm     , E_Int*   ipt_nidom_loc      , E_Int& it_bloc             , E_Int*   ipt_nisdom_residu,
                         E_Int*   ipt_it_lu_ssdom, E_Int*   ipt_it_target_ssdom, E_Int*   ipt_it_target_old , E_Int*   ipt_no_lu, E_Int*   param_int );

  void cprdu3s1_(E_Int& nd               , E_Int& nitcfg     , E_Int& nssiter            , E_Int& neq              ,  E_Int& ndimdx       , E_Int& ndim_rdm  , E_Int& nitrun, 
                 E_Int& iflw             , E_Int& iles       , E_Int& ithread            , E_Int& Nbre_thread_actif,  E_Int& omp_mode     ,  E_Int& mx_ssdom_lu  ,
                 E_Int& it_bloc          , E_Float& epsi     , E_Int*  size_ssdom        , E_Int* ipt_nisdom_residu,
                 E_Int*  nijk            , E_Int* ijkv       , E_Int*  ijkv_lu           , 
                 E_Int* ind_loop         ,
                 E_Int* ipt_it_lu_ssdom  , E_Int* ipt_it_target_ssdom,  E_Int* ipt_it_target_old, E_Int* ipt_it_temp_ssdom,  E_Int* ipt_no_lu, 
                 E_Float* iptdrodm       , E_Float* rdm  );


  void indice_boucle_lu_(E_Int& ndom , E_Int& ndsdm           , E_Int& nsdom_lu, E_Int& icp ,
                        E_Int* ind_dm, E_Int*  thread_topology, E_Int*  ind_sdm);

  void copy_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,   E_Float* iptro, E_Float* stock, E_Int& ind, E_Int& nzone);


  void copyrk3local_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,   E_Float* iptro, E_Float* stock, E_Int& ind, E_Int& pos, E_Int& taille);

  void copy_rk3local_(E_Int* ipt_param_int , E_Int* ind_loop, E_Float* iptro, E_Float* stock, E_Int& ind,E_Int& taille);

  void copy_rk3localpara_(E_Int* ipt_param_int , E_Int* ind_loop,E_Int* ind_loop_, E_Float* iptro, E_Float* stock, E_Int& ind,E_Int& taille);

  void extrap_rk3localpara_(E_Int* ipt_param_int , E_Int* ind_loop, E_Float* iptro);

  void copy_rk3locallist_(E_Int* ipt_param_int , E_Float* iptro, E_Float* stock, E_Int* listDnr, E_Int& nbpts);
 
  void verrou_c_(E_Int* lok , E_Int& type );

  void flush_real_(  E_Int& size, E_Float* tab);

  void interpolation_(E_Int& idir , E_Int* ipt_param_int , E_Float* ipt_param_real ,E_Int* ind_loop, E_Float* iptro_tmp , E_Float* iptro);

  void interpolation2_(E_Int& idir , E_Int* ipt_param_int , E_Float* ipt_param_real ,E_Int* ind_loop, E_Float* iptro_tmp , E_Float* iptro);

  void interpolation4_(E_Int& idir , E_Int* ipt_param_int , E_Float* ipt_param_real ,E_Int* ind_loop, E_Float* iptro_tmp , E_Float* iptro);

  void interpolation3_(E_Int& idir , E_Int* ipt_param_int , E_Float* ipt_param_real ,E_Int* ind_loop, E_Float* iptro_tmp , E_Float* iptro, E_Float* stock, E_Int& pos);

  void interprk3local_(E_Int& idir, E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Float* stock, E_Float* drodmstock, E_Int& nzone, E_Int& pos, E_Int& taille);

  void interp_rk3local_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Float* stock, E_Float* drodmstock,E_Int& taille);

  void interp_rk3local4_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Float* stock, E_Float* drodmstock,E_Int& taille, E_Float& coeff);

  void interp_rk3local4para_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Int* ind_loop_, E_Float* stock, E_Float* drodmstock,E_Float* iptro, E_Int& taille, E_Float& coeff);

  void interprk3local2_(E_Int& idir, E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Float* stock, E_Float* iptro, E_Int& nzone, E_Int& pos, E_Int& taille,E_Int& nstep);

  void interp_rk3local2_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Float* stock,E_Float* stock2, E_Float* iptro, E_Int& dir, E_Int& taille,E_Int& nstep);

  void interp_rk3local3_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop, E_Float* stock,E_Float* stock2, E_Float* iptro, E_Int& dir, E_Int& taille,E_Int& nstep);

  void interp_rk3local3para_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float*coe, E_Int* ind_loop,  E_Int* ind_loop_, E_Float* stock,E_Float* stock2, E_Float* iptro, E_Int& dir, E_Int& taille,E_Int& nstep,E_Int& NoD);

  void interp_rk3local3list_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Int*listRcv,E_Float* stock,E_Float* stock2, E_Float* iptro, E_Int& nbpts,E_Int& nstep);

  void conservativite32_(E_Int& idir , E_Int* ipt_param_int , E_Float* ipt_param_real, E_Int* ipt_param_int2, E_Int* ind_loop, E_Int& imax, E_Float* iptro,  E_Float*drodmzone, E_Float*drodm, E_Float* stock, E_Float*coe,E_Int& nstep, E_Int& nzone);

  void conservrk3local_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop, E_Float*drodm ,  E_Float*coe , E_Float* constk, E_Int& pos , E_Int& taille , E_Int& nstep, E_Int& ind);

  void conservrk3local2_(E_Int& idir , E_Int* ipt_param_int , E_Int* ipt_param_int2,E_Float*coe,E_Float*coe2, E_Float* ipt_param_real, E_Int* ind_loop, E_Float* iptro ,  E_Float* iptcstk1 , E_Float* iptcstk2, E_Int& taille , E_Int& nstep);

  void conservrk3local3_(E_Int* ipt_param_int , E_Int* ind_loop, E_Float*drodm ,  E_Float*coe , E_Float* constk, E_Int& taille , E_Int& nstep, E_Int& ind);

  void conservrk3local3para_(E_Int* ipt_param_int , E_Int* ind_loop, E_Int* ind_loop_, E_Float*drodm ,  E_Float*coe , E_Float* constk, E_Int& taille , E_Int& nstep, E_Int& ind);

  void conservrk3local4_(E_Int* ipt_param_int , E_Int* ipt_param_int2,E_Float*coe,E_Float*coe2, E_Float* ipt_param_real, E_Int* ind_loop, E_Float* iptro ,  E_Float* iptcstk1 , E_Float* iptcstk2, E_Int& taille , E_Int& nstep);

  void conservrk3local4para_(E_Int* ipt_param_int , E_Int* ipt_param_int2,E_Float*coe,E_Float*coe2, E_Float* ipt_param_real, E_Int* ind_loop,  E_Int* ind_loop_, E_Int* ind_loop__, E_Float* iptro ,  E_Float* iptcstk1 , E_Float* iptcstk2, E_Int& tailleR , E_Int& tailleD ,E_Int& nstep , E_Int& dir,E_Int& pt1,E_Int& pt2 ,E_Int& pt3,E_Int* transfo,E_Int& ratio1,E_Int& ratio2 ,E_Int& ratio3);

  void copyfluxrk3local_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,  E_Float*drodm, E_Float* stock, E_Int& ind, E_Int& pos,E_Int& taille);

  void copyflux_rk3local_(E_Int* ipt_param_int , E_Int* ind_loop, E_Float*drodm, E_Float* stock, E_Int& ind,E_Int& taille);

  void copyflux_rk3localpara_(E_Int* ipt_param_int , E_Int* ind_loop,E_Int* ind_loop_, E_Float*drodm, E_Float* stock, E_Int& ind,E_Int& taille);

  void copyflux_rk3locallist_(E_Int* ipt_param_int , E_Float*drodm, E_Float* stock, E_Int* listDnr, E_Int& npts, E_Int& ind);

  void copyflux_rk3local2_(E_Int* ipt_param_int , E_Int* ind_loop,  E_Float*drodm, E_Float* stock, E_Int& ind,E_Int& taille);

  void copyflux_rk3local2para_(E_Int* ipt_param_int , E_Int* ind_loop, E_Int* ind_loop_,  E_Float*drodm, E_Float* stock, E_Int& ind,E_Int& taille);

  void copyfluxrk3local2_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,  E_Float*drodm, E_Float* stock, E_Int& ind, E_Int& nzone,E_Int& taille);


  void constantinescu_(E_Float* iptro,E_Int* ind_loop,E_Float* stock,  E_Float*drodm, E_Float* coe,E_Int* ipt_param_int,E_Float* ipt_param_real,E_Int& nstep,E_Int& nzone);

  void constantinescurk3_(E_Int& idir,E_Int* ipt_param_int,E_Float* ipt_param_real,E_Int* ind_loop, E_Float* iptro,E_Float* iptstk,E_Int& nzone);

  void constantinescurk32_(E_Int& idir,E_Int* ipt_param_int,E_Float* ipt_param_real,E_Int* ind_loop, E_Float* iptro,E_Float* iptstk,E_Int& nzone);

   
  void copynuma_( E_Int*  ind_loop,  E_Int& ni,  E_Int& nj, E_Int& shift, E_Int& ific, E_Int& jfic, E_Int& kfic,
                 E_Float* target       , E_Float* src );

  void copyflux_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,   E_Float*drodm, E_Float* stock, E_Int& ind, E_Int& nzone);

  void remp_cellfictives_(E_Int* ipt_param_intDnr , E_Int* ipt_param_intRcv , E_Float*ropDnr, E_Float* ropRcv, E_Int* listDnr , E_Int* listRcv, E_Int& nbpts);

  void remp_cellfictivespara_(E_Int* ipt_param_intDnr , E_Int* ipt_param_intRcv , E_Float*ropDnr, E_Float* ropRcv, E_Int* listDnr , E_Int* listRcv, E_Int& nbpts,E_Int& pt_deb, E_Int& pt_fin);

  void copyfluxrk3local_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,  E_Float*drodm, E_Float* stock, E_Int& ind, E_Int& pos,E_Int& taille);

  void copyfluxrk3local2_(E_Int& idir , E_Int* ipt_param_int , E_Int* ind_loop,  E_Float*drodm, E_Float* stock, E_Int& ind, E_Int& nzone,E_Int& taille);

  void switchvectors_(E_Float* iptro, E_Float* iptrotmp, E_Int* ipt_param_int);

  void initdrodm_(E_Int* ipt_param_int , E_Int* ind_loop,   E_Float* drodm, E_Float* drodm2);

  void calcul_cfl_(E_Int& ndom, E_Int* ipt_param_int , E_Float* ipt_param_real,  E_Int* ind_loop, E_Float* ipt_cfl,   E_Float* iptro, E_Float* xmut, E_Float* venti, 
		   E_Float* ti, E_Float* tj, E_Float* tk, E_Float* vol, E_Float* ipt_cfl_, E_Int& isconv, E_Int& isvisc, E_Int& isSound);

  void init_ventijk_( E_Int& ndo  , E_Int& nidom  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   ,
                      E_Int* ipt_param_int        , E_Float* ipt_param_real , 
                      E_Int* ipt_ijkv_sdm       ,
                      E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_topo_thread  , E_Int* ipt_lok, E_Int* ipt_topo_omp, E_Int* ipt_inddm_omp,
                      E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                      E_Float* ipti_df          , E_Float* iptj_df     , E_Float* iptk_df        , 
                      E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                      E_Float* iptx             , E_Float* ipty         , E_Float* iptz    );

  void copy_ventijk_( E_Int& ndo                , E_Int& ithread      , E_Int* ipt_param_int, E_Float* ipt_param_real , 
                      E_Float* iptx             , E_Float* ipty       , E_Float* iptz,
                      E_Int* ipt_ind_dm         , E_Int* ipt_inddm_omp,
                      E_Float* iptventi         , E_Float* iptventj   , E_Float* iptventk   ,  E_Float* iptvent_vertex);

  void post_( E_Int& ndo  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& neq_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real  ,  E_Float& tke      , E_Float& enst     ,  E_Float& compteur  ,
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptgrad   );

  void post_grad_( E_Int& ndo,  E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& neq_grad, E_Int& order,
              E_Int* ipt_param_int      , E_Float* ipt_param_real  ,
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptgrad   );


  void post_q_( E_Int& ndo  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ,  E_Float* iptenst ,  E_Float* iptvort   );

  void post_qprime_( E_Int& ndo  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   ,E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            , E_Float* iptro_m1    ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ   );

  void post_qprime_1rot_( E_Int& ndo,   E_Int& Nbre_thread_actif , E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   ,
              E_Int& order              ,  E_Int& dim_grad          , E_Int& var1, E_Int& var2       ,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            , E_Float* iptro_m1    ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ , E_Float* iptrot );


  void post_enst_( E_Int& ndo  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   ,E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ,  E_Float* iptenst ,  E_Float* iptvort   );

  void post_q_enst_( E_Int& ndo  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  ,E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , E_Float* iptQ,  E_Float* iptenst ,  E_Float* iptvort   );


  void post_drodt_( E_Int& ndo  , E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& flag,  E_Int& dim_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real   , 
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_topo_sock    , E_Int* ipt_lok         ,
              E_Float* iptro            , E_Float*rop_n         , E_Float* rop_n1         , E_Float*drodt );


  void viles_( E_Int& ndo  ,  E_Int& Nbre_thread_actif, E_Int& ithread    , E_Int& Nbre_socket,  E_Int& socket    , E_Int& mx_synchro   , E_Int& neq_grad,
              E_Int* ipt_param_int      , E_Float* ipt_param_real  ,
              E_Int* ipt_ijkv_sdm       ,
              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_topo  , E_Int*  ipt_ind_dm_omp        , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
              E_Float* iptro            ,
              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol          ,  E_Float* iptmut, E_Float* iptdist, E_Float* iptrot   );

  void navier_stokes_struct_( E_Int& ndo    , E_Int& Nbre_thread_actif,
                              E_Int& ithread        , E_Int& ithread_io       , E_Int& omp_mode, E_Int& layer_mode, E_Int& Nbre_socket, E_Int& socket     , E_Int& mx_synchro   , 
                              E_Int& lssiter_verif  , E_Int& lexit_lu, E_Int& nptpsi           , E_Int& nitcfg  , E_Int& nssiter   , E_Int& nitrun    , E_Int& first_it   , E_Int& nb_pulse   , E_Int&   flagCellN  ,
                              E_Int* ipt_param_int  , E_Float* ipt_param_real ,
                              E_Float& temps        , E_Int* ipt_tot,   
                              E_Int* ipt_ijkv_sdm       ,
                              E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_ind_dm_omp  , E_Int* ipt_topo_thread  , E_Int* ipt_lok,  E_Int* ipt_topo_omp, E_Float* timer_omp,
                              E_Float* krylov           , E_Float& norm_kry,
                              E_Float* ipt_cfl        ,
                              E_Float* iptx             , E_Float* ipty        , E_Float* iptz           , E_Float* iptCellN       , E_Float* iptCellN_IBC  ,
                              E_Float* iptro            , E_Float* iptro_m1    , E_Float* iptrotmp       , E_Float* iptro_ssiter   ,
                              E_Float* iptmut           ,
                              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
                              E_Float* ipti_df          , E_Float* iptj_df     , E_Float* iptk_df        , E_Float* iptvol_df      , 
                              E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                              E_Float* iptwig           , E_Float* iptstat_wig , E_Float* iptrot         ,
                              E_Float* iptdrodm         , E_Float* iptcoe      , E_Float* iptdelta       ,  E_Float* iptro_res , E_Float* iptsrc);


  void invlu_(                E_Int& ndo      , E_Int& nitcfg      , E_Int& nitrun   , E_Int*  param_int , E_Float* param_real,
                              E_Int* ipt_sdm            , E_Int* ipt_sdm_thread, E_Int& mjrnewton        ,
                              E_Float* iptrotmp         , E_Float* iptro_ssiter,
                              E_Float* iptdrodmin       , E_Float* iptdrodmout ,
                              E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , 
                              E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
                              E_Float* iptcoe           , E_Float* iptssor     , E_Int& size_ssor);

  void dp_dw_vect_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Int* ind_loop,
		   E_Float* iptrop, E_Float* iptvectin, E_Float* iptvectout, E_Int& size);

  void bvbs_wall_inviscid_jacob_(E_Int* ipt_param_int, E_Int* ind_loop, E_Int& idir, E_Int& neq_mtr, E_Float* ipttijk, E_Float* iptvect);

  void bvbs_wall_inviscid_d_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                 E_Float* iptventi  , E_Float* iptijk  , E_Float* iptro, E_Float* iptrod);

  void pre_bc_(E_Int* param_int, E_Float& signe, E_Int*ind_loop, E_Float* krylov, E_Float* rop);

  void normalisation_vect_( E_Float& normL2, E_Int* param_int, E_Int* ind_loop, E_Float* krylov);

  void norm_vect_( E_Int* param_int, E_Int* ind_loop, E_Float* krylov, E_Float& norm);

  void id_vect_(E_Int* param_int, E_Int* ind_loop,  E_Float* drodmd,  E_Float* krylov_out, E_Float* krylov_in, E_Int& size);

  void navier_stokes_struct_d_( E_Int& ndo     , E_Int& nidom       , E_Int& Nbre_thread_actif, E_Int& ithread    ,
				E_Int& omp_mode, E_Int& layer_mode  , E_Int& Nbre_socket      , E_Int& socket     , E_Int& mx_synchro   , 
				E_Int& lssiter_verif  , E_Int& nptpsi            , E_Int& nitcfg  , E_Int& nitrun     ,
				E_Int& first_it   ,  E_Int& nb_pulse  , E_Int&   flagCellN  , E_Int& mjr_dt,
				E_Int* ipt_param_int  , E_Float* ipt_param_real  ,
				E_Float& temps        , E_Int* ipt_tot,   
				E_Int* ipt_ijkv_sdm       ,
				E_Int* ipt_ind_dm_int     , E_Int* ipt_ind_dm_sock, E_Int* ipt_inddm_omp   , E_Int* ipt_topo_thread  , E_Int* ipt_lok         ,
				E_Int* ipt_topo_omp,
				E_Float* ipt_cfl        ,
				E_Float* iptx             , E_Float* ipty         , E_Float* iptz          , E_Float* iptCellN       ,
				E_Float* iptro_ssiter     , E_Float* iptro_ssiterd, E_Float* krylov_in     ,
				E_Float* iptmut           , E_Float* iptmutd     ,
				E_Float* ipti             , E_Float* iptj        , E_Float* iptk           , E_Float* iptvol         , 
				E_Float* ipti_df          , E_Float* iptj_df     , E_Float* iptk_df        , E_Float* iptvol_df      , 
				E_Float* iptventi         , E_Float* iptventj    , E_Float* iptventk       ,  
				E_Float* iptwig           , E_Float* iptstat_wig , E_Float* iptrot         ,
				E_Float* iptdrodm         , E_Float* iptdrodmd   , E_Float* iptcoe         ,
				E_Float* iptdelta         , E_Float* iptro_res );

  void scal_prod_(E_Int* ipt_param_int, E_Int* ind_loop, E_Float* iptvect1,
		  E_Float* iptvect2, E_Float& value);

  void vect_rvect_(E_Int* ipt_param_int, E_Int* ind_loop, E_Float* iptvect1,
  		   E_Float* iptvect2, E_Float& value, E_Float& normL2);

  void prod_mat_vect_(E_Int* ipt_param_int, E_Int* ind_loop, E_Float* iptkrylov,
		      E_Float* iptvecty, E_Float* iptdrodm, E_Int& num_krylov);

  void mjr_prim_from_cons_(E_Int* ipt_param_int, E_Float* ipt_param_real, E_Float* visco, E_Float* sa_real, E_Int* ind_loop, 
		   E_Float* roptmp, E_Float* ropssiter, E_Float* drodm);
  


   void     bvbs_extrapolate_( E_Int& idir , E_Int& lrhs    , E_Int& eq_deb     ,E_Int* param_int , E_Int* ind_loop  , E_Float& nutildeinf,  E_Float* iptro);
   void     bvbs_extrapolate_d_( E_Int& idir , E_Int& lrhs    , E_Int& eq_deb     ,E_Int* param_int , E_Int* ind_loop  , E_Float& nutildeinf,  E_Float* iptro, E_Float* iptrod );

   void     bvbs_periodique_(   E_Int& idir  , E_Int& lrhs    , E_Int* param_int  , E_Int* ind_loop, E_Float* iptro);
   void     bvbs_periodique_d_( E_Int& idir  , E_Int& lrhs    , E_Int* param_int  , E_Int* ind_loop, E_Float* iptro, E_Float* iptrod);

   void     bvbs_periodique_azimuthal_(   E_Int& idir  , E_Int& lrhs    , E_Int* param_int  , E_Int* ind_loop, E_Float* iptro, E_Float* iptdata);
   void     bvbs_periodique_azimuthal_d_( E_Int& idir  , E_Int& lrhs    , E_Int* param_int  , E_Int* ind_loop, E_Float* iptro,  E_Float* iptrod, E_Float* iptdata);

   void     bvbs_wall_inviscid_(   E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                   E_Float* iptventi  , E_Float* iptijk  , E_Float* iptro);
   void     bvbs_wall_inviscid_d_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                   E_Float* iptventi  , E_Float* iptijk  , E_Float* iptro, E_Float* iptrod);

   void     bvbs_wallmodel_(   E_Int& idir      , E_Int& lrhs      , E_Int& neq_mtr, E_Float& mobile_coef, E_Float& c4   , E_Float& c5, E_Float& c6, 
                               E_Int* param_int , E_Int* ind_loop  ,
                               E_Float* iptventi, E_Float* iptijk  , E_Float* iptro, E_Float* xmut);

   void     bvbs_wallexchange_(   E_Int& idir      , E_Int& lrhs      , E_Int& neq_mtr, E_Float& mobile_coef, E_Float& c4   , E_Float& c5, E_Float& c6, 
                               E_Int* param_int , E_Float* param_real, E_Int* ind_loop  ,
                               E_Int& sizetab   , E_Int* inc_bc,
                               E_Float* iptx    , E_Float* ipty    , E_Float* iptz ,
                               E_Float* iptventi, E_Float* iptijk  , E_Float* iptro, E_Float* xmut, E_Float* moy, E_Float* wmles);

   void     bvbs_wall_viscous_adia_(   E_Int& idir      , E_Int& lrhs      , E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                       E_Float* iptventi, E_Float* iptijk  , E_Float* iptro);
   void     bvbs_wall_viscous_adia_d_( E_Int& idir      , E_Int& lrhs      , E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                       E_Float* iptventi, E_Float* iptijk  , E_Float* iptro, E_Float* iptrod);

   void     bvbs_wall_viscous_isothermal_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Float& mobile_coef, E_Int* param_int ,E_Int* ind_loop  ,
                                              E_Float* param_real, E_Float* iptventi, E_Float* iptijk, E_Float* iptro,
                                              E_Float* state1    , E_Int& size_data , E_Int* inc_bc  , E_Int& size_work);

   void     bvbs_wall_viscous_transition_(   E_Int& idir      , E_Int& lrhs      ,  E_Int& nstep,   E_Int& neq_mtr   , E_Float& mobile_coef, E_Float* random, 
					     E_Int& size_data , E_Int* inc_bc    ,  E_Int* param_int , E_Int* ind_loop     , E_Float* param_real,
                                             E_Float* iptx    , E_Float* ipty     , E_Float* iptz  ,
                                             E_Float* iptventi, E_Float* iptijk   , E_Float* iptro);
 
   void     bvbs_wall_viscous_transition_d_( E_Int& idir     , E_Int& lrhs  ,  E_Int& nstep    , E_Int& neq_mtr  , E_Float& mobile_coef, E_Float* random, 
					     E_Int& size_data, E_Int* inc_bc,  E_Int* param_int, E_Int* ind_loop     , E_Float* param_real,
                                             E_Float* iptx    , E_Float* ipty     , E_Float* iptz ,
                                             E_Float* iptventi, E_Float* iptijk   , E_Float* iptro,  E_Float* iptrod);


   void     bvbs_inflow_supersonic_(   E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                       E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                       E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);
   void     bvbs_inflow_supersonic_d_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                       E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                       E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod, E_Float* state);

   void     bvbs_farfield_(   E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                              E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                              E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);
   void     bvbs_farfield_d_( E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                              E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                              E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod, E_Float* state);

   void     bvbs_outflow_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                              E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                              E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);
   void     bvbs_outflow_d_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                              E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                              E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod, E_Float* state);

   void     bvbs_inflow_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                             E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                             E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* state);
   void     bvbs_inflow_d_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                             E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                             E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod, E_Float* state);


   void     bvbs_outpres_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr,  E_Int* param_int ,E_Int* ind_loop  ,
                              E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                              E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* data_pres,
			      E_Int& size_data, E_Int* inc_bc);
   void     bvbs_outpres_d_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr,  E_Int* param_int ,E_Int* ind_loop  ,
                              E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                              E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod, E_Float* data_pres,
			      E_Int& size_data, E_Int* inc_bc);

//   void     bvbs_updatepressure_(  E_Int& idir  , E_Int& ithread, E_Int* ind_avg, E_Int* ind_mjr, E_Int* param_in ,E_Int* ind_loop ,
//                            E_Float* param_real, E_Float* ipty   , E_Float* iptz  , E_Float* iptro, E_Float* data_pres,
//			    E_Int& size_data, E_Float* vteta, E_Float* roteta, E_Int& inc_bc);

   void     bvbs_inflow_newton_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                    E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                    E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, 
                                    E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	            E_Int& size_data   , E_Int* inc_bc     , E_Int& size_work);
   void     bvbs_inflow_newton_d_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                    E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                    E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod,
                                    E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	            E_Int& size_data   , E_Int* inc_bc     , E_Int& size_work);

   void     bvbs_inflow_lund_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                  E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                  E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, 
                                  E_Float* roIn      , E_Float* roLund   , E_Float* param_lund, 
	         	          E_Int& size_data   , E_Int* inc_bc     , E_Int& size_work);


   void     bvbs_inflow_fich_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                  E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                  E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, 
                                  E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	          E_Int& size_data   , E_Int* inc_bc     , E_Int& size_work);
   void     bvbs_inflow_fich_d_(  E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                  E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                  E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, E_Float* iptrod,
                                  E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	          E_Int& size_data   , E_Int* inc_bc     , E_Int& size_work);

   void     bvbs_inflow_supersonic_fich_(   E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                                       E_Float* param_real, E_Float& c4   , E_Float& c5, E_Float& c6,
                                       E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro, 
                                       E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6, 
	         	               E_Int& size_data   , E_Int* inc_bc     , E_Int& size_work);

   void     bvbs_injmfr_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                            E_Float* param_real, E_Float& c4   ,     E_Float& c5, E_Float& c6,
                            E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro,
                            E_Float* state1    , E_Float* state2   , E_Float* state3, E_Float* state4, E_Float* state5, E_Float* state6,
                            E_Int& size_data   , E_Int* inc_bc  );

   void    bvbs_outmfr_(    E_Int& idir        , E_Int& lrhs      ,  E_Int& neq_mtr, E_Int* param_int ,E_Int* ind_loop  ,
                           E_Float* param_real, E_Float& c4   ,     E_Float& c5, E_Float& c6,
                           E_Float* iptventi  , E_Float* iptijk   , E_Float* iptro,
                           E_Float* state1    ,
                           E_Int& size_data   , E_Int* inc_bc  );

   void indice_cl_sdm_( E_Int& dir    , E_Int& npass  , E_Int& lskip , E_Int& typbc,
                        E_Int& ific   , E_Int& kfic   , 
                        E_Int* ijkv   , E_Int* ind_fen, E_Int* ind_dm, E_Int* ind_dm_thread, 
                        E_Int* ind_CL , E_Int* ind_CL119);

   void     correct_coins_( E_Int& ndo ,  E_Int* param_int, E_Int* ind_loop, E_Float* iptro);

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


   void     move_domx_( E_Int& ndo ,  E_Int* param_int, E_Float* param_real, E_Int* ind_loop, 
                        E_Float* iptx  , E_Float* ipty , E_Float* iptz, 
                        E_Float* iptx0 , E_Float* ipty0, E_Float* iptz0  );

   void     mjr_ale_3dhomocart_( E_Int& ndo  , E_Int* param_int   , E_Float* param_real,
                        E_Int& socket        , E_Int& Nbre_socket , E_Int& thread_sock , 
                        E_Int& thread_parsock, E_Int* ind_socket  , E_Int* topo,
                        E_Float* iptx  , E_Float* ipty , E_Float* iptz, 
                        E_Float* ipti  , E_Float* iptj , E_Float* iptk ,
                        E_Float* ipti0 , E_Float* iptj0, E_Float* iptk0, E_Float* iptvol ,E_Float* venti, E_Float* ventj, E_Float* ventk);


  void     dpssiter_(  E_Int& nitrun , E_Int& neq , E_Int& nssiter, E_Int& iflw, E_Int& iles, E_Int& lft, E_Int& iverb,  char*, E_Int& size_name, E_Float* rdm, E_Float* cvg_ptr);

  void     conv2pytree_(E_Float* cvg_pt, E_Int& nitrun, E_Int& neq, E_Int* LastRec, char* name, E_Int& size_name, E_Int& lft, E_Int& nrec, E_Int& nd, E_Int* Itnum, E_Float* Res_L2, E_Float* Res_oo, E_Float* Res_L2_diff, E_Float* Res_oo_diff);

  void bceffort_( E_Int& ndo        ,  E_Int& ithread     ,  E_Int* param_int, E_Float* param_real,  E_Int* param_int_eff,
                  E_Int* ind_loop   ,  E_Float* effort_omp,  E_Float* xyz_ref,
                  E_Float*  iptro   ,  E_Float* iptflu    ,  E_Float* iptwig , E_Float* iptmut,   
                  E_Float* iptx     ,  E_Float* ipty      ,  E_Float* iptz   ,
                  E_Float* ipti     ,  E_Float* iptj      ,  E_Float* iptk   , E_Float* iptvol,
                  E_Float* venti    ,  E_Float* ventj     ,  E_Float* ventk );

  void sfd_( E_Int* param_int, E_Float* param_real, E_Int& nitrun, E_Int* ind_loop, E_Float* ipt_CL, E_Float* iptrof, E_Float* iptcoe, E_Float* iptvol );

  void cp_conservatif_( E_Int* param_int, E_Float* param_real, E_Int& nthreadmax, E_Int* ind_loop, E_Float* ipt_CL, E_Float* iptvol, E_Float* iptcellN, E_Float* debit  );

  void corr_conservatif_( E_Int* param_int, E_Float* param_real, E_Int& nthreadmax, E_Int* ind_loop, E_Float* ipt_CL, E_Float* iptvol, E_Float* iptcellN, E_Float* debit  );

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

 void mj_lund_planrecyl_(E_Int& nd, E_Int& idir, E_Int* param_int , E_Int* ind_fen,  E_Int* ind_dm_th, E_Int* inc_bc, E_Int& size_data, E_Float* lund_param, E_Float* rop, E_Float* rof2 );

 void mj_wallmodel_plan_(E_Int& nd, E_Int& idir, E_Int* param_int , E_Int* ind_fen,  E_Int* ind_dm_th, E_Int* inc_bc, E_Int& size_data, E_Int& lprint, E_Float* lund_param, E_Float* rop, E_Float* moy,  E_Float* snap  );

 void cp_debit_ibm_(E_Int& nd, E_Int& idir, E_Int& neq_mtr, E_Int& ithread, E_Int& Nthread_max,  E_Int& nitcfg, E_Int* param_int , E_Float* param_real, 
                    E_Int& size_fen, E_Int* facelist,  E_Float* rop, E_Float* iptijk, E_Float* iptvol ,    E_Float* flux   );

 void cp_corr_debit_ibm_(E_Int& nd, E_Int& idir, E_Int& neq_mtr, E_Int& ithread, E_Int& Nthread_max,  E_Int& nitcfg, E_Int* param_int , E_Float* param_real, 
                    E_Int& size_fen, E_Int* facelist,  E_Float* rop, E_Float* iptijk, E_Float* coe ,  E_Float* drodm ,  E_Float* flux   );

 void corr_debit_ibm_(E_Int& nd, E_Int& idir, E_Int& neq_mtr, E_Int& ithread, E_Int& Nthread_max, E_Int* param_int , E_Float* param_real, E_Int& ipass, E_Int& fam, E_Int& nstep, 
                     E_Float* amor, E_Int& size_fen, E_Int* facelist,  E_Float* rop, E_Float* iptijk, E_Float* iptvol ,  E_Float* celln,  E_Float* flux   );

 void corr_bilan_ibm_(E_Int& nd, E_Int& idir, E_Int& neq_mtr, E_Int& ithread, E_Int& Nthread_max, E_Int& nitcfg, E_Int* param_int , E_Float* param_real,
                      E_Int& size_fen, E_Int* facelist,  E_Float* drodm, E_Float* iptijk, E_Float* iptcoe , E_Float* flux   );

  }
