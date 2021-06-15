E_Int nd_subzone = 0;
for (E_Int nd = 0; nd < nidom; nd++){
  //mise a jour Nombre sous_iter pour implicit (simplification gestion OMP)
  if( (param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1 || lssiter_verif ==1 ) && param_int[nd][ ITYPZONE] != 4){
    E_Int* ipt_nisdom_residu   =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter;                //nisdom_residu(nssiter)
    E_Int* ipt_it_bloc         =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter*2;    //it_bloc(nidom)

    if(ipt_nisdom_residu[nitcfg-1] != 0) ipt_it_bloc[0] +=1;
  }

  //E_Int* ipt_nidom_loc= ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;


  if( (param_int[ nd ][ NB_RELAX ] > 1 || param_int[ nd ][ LU_MATCH ]==1) && param_int[nd][ ITYPZONE] != 4){
#ifdef _OPENMP
    E_Int  Nbre_thread_actif = __NUMTHREADS__;
#else
    E_Int Nbre_thread_actif = 1;
#endif
    //E_Int nfic_ij = param_int[ nd ][ NIJK + 3 ];
    //E_Int nfic_k  = param_int[ nd ][ NIJK + 4 ];
    E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
    E_Int  nb_subzone    = ipt_nidom_loc [nitcfg-1];                                           //nbre sous-zone a la sousiter courante

    for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++){
      if (omp_mode == 1){
	E_Int       Ptomp = param_int[nd][PT_OMP];
	E_Int  PtrIterOmp = param_int[nd][Ptomp];   
	E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

	for (E_Int i = 0; i < Nbre_thread_actif; i++){
	  if (param_int[nd][PtZoneomp + i] != - 2){
	    E_Int* ipt_ind_dm_thread = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (param_int[nd][PtZoneomp + i]) * 6;
	    E_Int indice             = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
	    ipt_ssor_size[indice]=(ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] +1 +2*param_int[nd][ NIJK +3 ])*
	      (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] +1 +2*param_int[nd][ NIJK +3 ])*
	      (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] +1 +2*param_int[nd][ NIJK +4 ]);
	  }
	  else { E_Int indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
	    ipt_ssor_size[indice]=0;
	  }
	}
      }
      else { //ompmode = 0
	E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
	E_Int ipt_topology_socket[3];
	E_Int ipt_ind_dm_thread[6];
	E_Int lmin = 10;
	if (param_int[nd][ITYPCP] == 2) lmin = 4;

	for (E_Int i = 0; i < Nbre_thread_actif; i++){
	  E_Int ithread = i + 1;
	  E_Int indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;

	  indice_boucle_lu_(nd, ithread, Nbre_thread_actif, lmin,
			    ipt_ind_dm_loc,
			    ipt_topology_socket, ipt_ind_dm_thread);

	  if (ithread > ipt_topology_socket[0]*ipt_topology_socket[1]*ipt_topology_socket[2]){
	    ipt_ssor_size[indice]=0;
	    //break;
	  }
	  else{
	    ipt_ssor_size[indice]=(ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] +1 +2*param_int[nd][ NIJK +3 ])*
	      (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] +1 +2*param_int[nd][ NIJK +3 ])*
	      (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] +1 +2*param_int[nd][ NIJK +4 ]);
	  }
	}
      }// ompmode
    }//loop subzone
  }  // if relax
 } // loop zone
