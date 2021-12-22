   E_Int* ipt_ind_dm_thread; 
   E_Int nd_current =0;


   for (E_Int nd = 0; nd < nidom; nd++)
     {
       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
       E_Int  nb_subzone    = ipt_nidom_loc [nitcfg-1];                                           //nbre sous-zone a la sousiter courante

       E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
       E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
       E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
       E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;

       for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
	 {
	   E_Int ndo   = nd;

	   E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

	   E_Int* ipt_ind_dm_thread;
	   if (omp_mode == 1)
	     { 
	       E_Int       Ptomp = param_int[nd][PT_OMP];
	       E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
	       E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

	       Nbre_thread_actif_loc = param_int[nd][ PtZoneomp  + Nbre_thread_actif ];
	       ithread_loc           = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
	       ipt_ind_dm_thread     = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;

	       if (ithread_loc == -1) { continue;}
	     }
	   else
	     { 
	       E_Int* ipt_topology_socket = ipt_topology       + (ithread-1)*3;
	       E_Int* ipt_ind_dm_socket   = ipt_ind_dm_omp     + (ithread-1)*12;
	       ipt_ind_dm_thread   = ipt_ind_dm_socket  +6;


	       E_Int lmin = 10;
	       if (param_int[nd][ITYPCP] == 2) lmin = 4;

	       indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, lmin,
				 ipt_ind_dm_loc,
				 ipt_topology_socket, ipt_ind_dm_thread);
	     }
	   if (autorisation_bc[nd] == 1)
	    {
	       E_Int ierr = BCzone(nd, lrhs , nitcfg_stk, lcorner, param_int[nd], param_real[nd], npass,
				   ipt_ind_dm_loc, ipt_ind_dm_thread, 
				   ipt_ind_CL_thread, ipt_ind_CL119,  ipt_ind_CLgmres, ipt_shift_lu,
				   iptro_CL[nd] , ipti[nd]            , iptj[nd]    , iptk[nd]       ,
				   iptx[nd]     , ipty[nd]            , iptz[nd]    ,
				   iptventi[nd] , iptventj[nd]        , iptventk[nd], iptro_CL[nd], iptmut[nd]);

               //if(nd>=274){  a ne pas appliquer si zone IBM classique:  overlap coin OK
               if(nd>=0){
	       correct_coins_(nd,  param_int[nd], ipt_ind_dm_thread , iptro_CL[nd]);
               }

	    }//autorisation


           //Reinitialisation verrou omp
           //
           E_Int l =  nd_current*mx_synchro*Nbre_thread_actif  + (ithread_loc-1)*mx_synchro;
           for (E_Int i = 0;  i < mx_synchro ; i++) { ipt_lok[ l + i ]  = 0; }
           nd_current +=1;

	 }//loop souszone
     }//loop zone
