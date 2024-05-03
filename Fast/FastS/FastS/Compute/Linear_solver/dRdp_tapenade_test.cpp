ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;
nb_subzone    = ipt_nidom_loc [nitcfg-1];

for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
  {
    E_Int ndo   = nd;
    E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp;

    if (omp_mode == 1)
      { 
	E_Int mx_sszone = mx_nidom/nidom;
	E_Int shift_omp = param_int[nd][ PT_OMP ] + mx_sszone*nd_subzone*(Nbre_thread_actif*7+4);

	Nbre_thread_actif_loc = param_int[nd][ shift_omp  + Nbre_thread_actif ];
	ithread_loc           = param_int[nd][ shift_omp  +  ithread -1       ] +1 ;
	ipt_topo_omp          = param_int[nd] + shift_omp +  Nbre_thread_actif + 1;
	ipt_inddm_omp         = param_int[nd] + shift_omp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;

	if (ithread_loc == -1) {continue;}
      }

    E_Int* ipt_lok_thread;
    ipt_lok_thread   = ipt_lok   + nd_current*mx_synchro*Nbre_thread_actif;

    E_Int* ipt_ind_dm_loc         = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
    E_Float* ipt_cfl_thread       = ipt_cfl         + (ithread-1)*3+ ndo*3*Nbre_thread_actif;
    E_Float* iptCellN_loc; E_Int flagCellN;
    if (iptCellN[nd] == NULL) { flagCellN = 0; iptCellN_loc = iptro[nd];}
    else                      { flagCellN = 1; iptCellN_loc = iptCellN[nd]; }

    indice_boucle_lu_(ndo, ithread_loc, Nbre_thread_actif_loc, lmin,
		      ipt_ind_dm_loc,
		      ipt_topology_socket, ipt_ind_dm_omp_thread );

    indice_boucle_lu_(ndo, socket , Nbre_socket, lmin,
		      ipt_ind_dm_loc,
		      ipt_topology_socket, ipt_ind_dm_socket );

    navier_stokes_struct_d_( ndo, nidom, Nbre_thread_actif_loc, ithread_loc, omp_mode, layer_mode, Nbre_socket, socket, mx_synchro , 
                             lssiter_verif, nptpsi, nitcfg, nitrun, first_it, nb_pulse, flagCellN,
			     param_int[nd] , param_real[nd] ,
			     temps               , ipt_tot       ,
			     ipt_ijkv_sdm_thread , ipt_ind_dm_loc, ipt_ind_dm_socket, ipt_ind_dm_omp_thread, ipt_topology_socket, ipt_lok_thread, ipt_topo_omp, ipt_inddm_omp,
			     ipt_cfl_thread      ,
			     iptx[nd]                , ipty[nd]                , iptz[nd]            , iptCellN_loc     ,
			     iptro_ssiter[nd]        , rop_ssiter_d            , krylov_in           ,
			     iptmut[nd]              ,  ipt_mutd  + shift_mu   , 
			     ipti[nd]                , iptj[nd]                , iptk[nd]            , iptvol[nd]       ,
			     ipti0[nd]               , iptj0[nd]               , iptk0[nd]           , iptvol_df[nd]    ,
			     iptventi[nd]            , iptventj[nd]            , iptventk[nd]        ,
			     iptwig   + shift_wig    , iptstat_wig + shift_wig , iptrot+ shift_wig   ,
			     krylov_in   , ipt_drodmd + shift_zone , iptcoe  + shift_coe     ,
			     iptdelta[nd]            , iptro_res[nd]     );

    nd_current++;

  }
