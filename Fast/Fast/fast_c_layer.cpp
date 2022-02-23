if(lexit_lu ==0 && layer_mode==1)
  {   

    //Swap (call to setInterpTransfer)

    //E_Float trans_begin = omp_get_wtime();


    E_Int cycl;
    E_Float deb_dtlocal;
    E_Float tmps;
    
    if (param_int[0][EXPLOC] == 0)
      {
	K_FASTC::setInterpTransfersFast(iptro_CL, vartype, param_int_tc, param_real_tc , param_int, param_real,
					linelets_int, linelets_real, it_target, nidom, ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage);
      }//test exploc==0

 
    if (param_int[0][EXPLOC] == 1 ) 
      {
	K_FASTS::dtlocal2para_c(iptro, iptrotmp, param_int_tc, param_real_tc, param_int, param_real, iptdrodm, iptcoe, stock, drodmstock, constk, nitcfg, omp_mode, taille_tabs, nidom);

	for (E_Int nd=0; nd < nidom; nd++)
	  {
	    cycl = nssiter/param_int[nd][LEVEL];
	    if (nitcfg%cycl == 0 and nitcfg != nssiter)
	      {
		E_Float* ptsave  = iptro[nd]; 
		iptro[nd] = iptrotmp[nd]; 
		iptrotmp[nd] = ptsave;
	      } 
	  }

	K_FASTC::setInterpTransfersFast(iptro_CL, vartype, param_int_tc, param_real_tc , param_int, param_real,
					linelets_int, linelets_real, it_target   , nidom        , ipt_timecount, mpi       , nitcfg, nssiter, rk, exploc, numpassage);


	// achtung
	//     recup3para_c(iptro, param_int_tc, param_real_tc, param_int, stock, nitcfg, omp_mode, taille_tabs, nidom);

	numpassage=2;
	if (nitcfg%2==0)
	  {
	    K_FASTC::setInterpTransfersFast(iptro_CL, vartype, param_int_tc, param_real_tc , param_int, param_real,
					    linelets_int, linelets_real, it_target   , nidom        , ipt_timecount, mpi       , nitcfg, nssiter, rk, exploc, numpassage);
	  }
 

      } // fin boucle test dtlocal

  


    //
    //E_Float trans_end = omp_get_wtime();
    //E_Float trans_duree = trans_end - trans_begin;
    //cout << "tps_trans  : "<< trans_duree  << endl;
    //
    //
    //Apply BC (parcour Zones) + reinitialisation verrou pour calcul rhs
    //
    //
    
    E_Int autorisation_bc[nidom];
    E_Int nitcfg_stk = nitcfg;
 
    for (E_Int nd = 0; nd < nidom; nd++)
      {
	// flag pour transfert optimiser explicit local
	autorisation_bc[nd]=0;
	if (exploc == 1) 
	  {
	    cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];
	  
	    if      (nitcfg_stk%cycl == cycl/2 -1       and cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}
	    else if (nitcfg_stk%cycl == cycl/2 + cycl/4 and cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}	    
	    else if (nitcfg_stk%cycl == cycl-1          and cycl != 4 ){ nitcfg = 1; autorisation_bc[nd] = 1;}
	    else if((nitcfg_stk%cycl == 1 or nitcfg_stk%cycl == cycl/2  or nitcfg_stk%cycl== cycl-1) and cycl == 4 ) { nitcfg = 1; autorisation_bc[nd] = 1; }
	  } 
	else {autorisation_bc[nd] = 1;}

	if(nitcfg==nitcfg_last-1 and param_int[nd][ITYPZONE] !=4  and param_int[nd][IFLOW] !=4)
	  {
	    ///mise a jour moyenne plan Lund si necessaire
	    E_Int pt_bcs = param_int[nd][PT_BC];
	    E_Int nb_bc  = param_int[nd][ pt_bcs ];
	    E_Float* ipt_data=NULL;
             E_Int nd_subzone =0;
	    for ( E_Int ndf = 0; ndf < nb_bc; ndf++ )
	      {
		E_Int pt_bc  = param_int[nd][pt_bcs+ 1 + ndf];

		E_Int idir   = param_int[nd][pt_bc + BC_IDIR];
		E_Int nbdata = param_int[nd][pt_bc + BC_NBDATA];
		E_Int bc_type= param_int[nd][pt_bc + BC_TYPE];

		if(bc_type==19) 
		  { 
		    E_Int* iptsize_data = param_int[nd] + pt_bc + BC_NBDATA + 1;
		    E_Int* ind_fen      = param_int[nd] + pt_bc + BC_FEN;
		    E_Int  inc_bc[3];
		    if ( idir <= 2 ) 
		      { inc_bc[0] = ind_fen[3] - ind_fen[2] + 1; // nombre element de la fenetre dans la direction J
			inc_bc[1] = ind_fen[2]; // debut indice j
			inc_bc[2] = ind_fen[4]; // debut indice k
		      }  
		    else if ( idir <= 4 )
		      { inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
			ind_fen[0]; // debut indice i
			inc_bc[2] = ind_fen[4]; // debut indice k
		      }  
		    else
		      {inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
			inc_bc[1] = ind_fen[0]; // debut indice i
			inc_bc[2] = ind_fen[2]; // debut indice j
		      }

		    if ( nbdata != 0 ) ipt_data = param_real[nd] + param_int[nd][pt_bcs + 1 + ndf + nb_bc];

		    E_Float* iptAvgPlanLund = ipt_data + 5*iptsize_data[0];
		    E_Float* iptParamLund   = ipt_data + 10*iptsize_data[0];

		    E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

		    mj_lund_planrecyl_(nd, idir, param_int[nd] , ind_fen, ipt_ind_dm_loc, inc_bc, iptsize_data[0],
				       iptParamLund,  iptro[nd], iptAvgPlanLund); 
		  }//lund
	      }//ndf
	  }//nit_last
      }//loop zone
 
    E_Int ndim=0;
    //return ndim;

    E_Int lrhs=0; E_Int lcorner=0; 
#pragma omp parallel default(shared)
    {
#ifdef _OPENMP
      E_Int  ithread           = omp_get_thread_num() +1;
      E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
      E_Int ithread = 1;
      E_Int Nbre_thread_actif = 1;
#endif

      //E_Int Nbre_socket   = NBR_SOCKET;             
      E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
      if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

      E_Int Nbre_thread_actif_loc, ithread_loc;
      if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
      else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}


      E_Int* ipt_ind_dm_thread; E_Int* ipt_topology_socket; E_Int* ipt_ind_dm_socket;
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
		  ipt_topology_socket = ipt_topology       + (ithread-1)*3;
		  ipt_ind_dm_socket   = ipt_ind_dm_omp     + (ithread-1)*12;
		  ipt_ind_dm_thread   = ipt_ind_dm_socket  +6;

		  E_Int lmin = 10;
		  if (param_int[nd][ITYPCP] == 2) lmin = 4;

		  indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, lmin,
				    ipt_ind_dm_loc,
				    ipt_topology_socket, ipt_ind_dm_thread);
		}

	      if (autorisation_bc[nd] == 1)
		{
		  if (param_int[nd][ITYPZONE ] !=4) // zone structuree
		    {
		      if(param_int[nd][IFLOW] != 4)
			{
			  E_Int ierr = K_FASTS::BCzone(nd, lrhs , nitcfg_stk, lcorner, param_int[nd], param_real[nd], npass,
						       ipt_ind_dm_loc, ipt_ind_dm_thread, 
						       ipt_ind_CL_thread, ipt_ind_CL119,  ipt_ind_CLgmres, ipt_shift_lu,
						       iptro_CL[nd] , ipti[nd]            , iptj[nd]    , iptk[nd]       ,
						       iptx[nd]     , ipty[nd]            , iptz[nd]    ,
						       iptventi[nd] , iptventj[nd]        , iptventk[nd], iptro_CL[nd], iptmut[nd]);
			  correct_coins_(nd,  param_int[nd], ipt_ind_dm_thread , iptro_CL[nd]);
			}
		      else
			{
			  //                     E_Int ierr = K_FASTLBM::BCzone(nd, lrhs , nitcfg_stk, lcorner, param_int[nd], param_real[nd], npass,
			  //                                                    ipt_ind_dm_loc, ipt_ind_dm_thread, 
			  //                                                    ipt_ind_CL_thread, ipt_ind_CL119,  ipt_ind_CLgmres, ipt_shift_lu,
			  //                                                    iptro_CL[nd] , ipti[nd]            , iptj[nd]    , iptk[nd]       ,
			  //                                                    iptx[nd]     , ipty[nd]            , iptz[nd]    ,
			  //                                                    iptventi[nd] , iptventj[nd]        , iptventk[nd], iptro_CL[nd]); 
			}
		    }
		  else
		    {
		      E_Int lrhs=0; E_Int lcorner=0;
		      K_FASTP::BCzone( nd, ithread, lrhs, lcorner,
				       param_int[nd], param_real[nd],
				       npass, temps,
				       ipt_ind_CL119  , 
				       ipt_ng_pe[nd]          , iptro_CL[nd]   ,
				       ipti[nd] , iptventi[nd], iptx[nd] , ipty[nd] , iptz[nd] );
		    }
		}//autorisation

	      //Reinitialisation verrou omp
	      //
	      E_Int l =  nd_current*mx_synchro*Nbre_thread_actif  + (ithread_loc-1)*mx_synchro;
	      for (E_Int i = 0;  i < mx_synchro ; i++) { ipt_lok[ l + i ]  = 0; }
	      nd_current +=1;

	    }//loop souszone
	}//loop zone

    }//fin zone omp


    nitcfg = nitcfg_stk;

    //     cout << "nidom= "<< nidom << endl;	 

    //cout << "coucou" << endl;



  }//test exit_lu

