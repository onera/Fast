    E_Int ndo   = nd;

    E_Int* ipt_topo_omp;

    ipt_topo_omp          = ipt_omp + pttask + 3 + Nbre_thread_actif ;

    if (ithread_loc == -1) {continue;}

    E_Float* ipt_cfl_thread= ipt_cfl         + (ithread-1)*3+ ndo*3*Nbre_thread_actif;

    E_Int* ipt_lok_thread;
    ipt_lok_thread   = ipt_lok   + ntask*mx_synchro*Nbre_thread_actif;

    E_Float* iptCellN_loc; E_Int flagCellN;
    if (iptCellN[nd] == NULL) { flagCellN = 0; iptCellN_loc = iptro[nd];}
    else                      { flagCellN = 1; iptCellN_loc = iptCellN[nd]; }

    indice_boucle_lu_(ndo, socket , Nbre_socket, lmin,
		      ipt_ind_dm_loc,
		      ipt_topology_socket, ipt_ind_dm_socket );

    navier_stokes_struct_d_( ndo, nidom, Nbre_thread_actif_loc, ithread_loc, omp_mode, layer_mode, Nbre_socket, socket, mx_synchro , 
                             lssiter_verif, nptpsi, nitcfg, nitrun, first_it, nb_pulse, flagCellN, mjr_dt,
			     param_int[nd] , param_real[nd] ,
			     temps               ,
			     ipt_ijkv_sdm_thread , ipt_ind_dm_loc, ipt_ind_dm_socket, ipt_ind_dm_thread, ipt_topology_socket, ipt_lok_thread, ipt_topo_omp,
			     ipt_cfl_thread      ,
			     iptx[nd]                , ipty[nd]                , iptz[nd]            , iptCellN_loc     ,
			     iptro_ssiter[nd]        , rop_ssiter_d            , krylov_in           ,
			     iptmut[nd]              , ipt_mutd + shift_mu,
			     ipti[nd]                , iptj[nd]                , iptk[nd]            , iptvol[nd]       ,
			     ipti0[nd]               , iptj0[nd]               , iptk0[nd]           , iptvol_df[nd]    ,
			     iptventi[nd]            , iptventj[nd]            , iptventk[nd]        ,
			     iptwig   + shift_wig    , iptstat_wig + shift_wig , iptrot+ shift_wig   ,
			     iptdrodm + shift_zone   , ipt_drodmd + shift_zone , iptcoe  + shift_coe     ,
			     iptdelta[nd]            , iptro_res[nd]     );
