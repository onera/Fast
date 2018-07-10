if (param_int[ nd ][ NB_RELAX ] > 1)
  {
    E_Int nfic_ij = param_int[ nd ][ NIJK + 3 ];
    E_Int nfic_k  = param_int[ nd ][ NIJK + 4 ];
    ipt_ssor_shift = iptssor[nd];

    ipt_ssor_size[ ithread - 1] = (ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] + 1 + 2 * nfic_ij) *
      (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] + 1 + 2 * nfic_ij) *
      (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] + 1 + 2 * nfic_k);

#pragma omp barrier

    if (omp_mode == 0)
      for (E_Int i = 0; i < ithread - 1; i++)
	ipt_ssor_shift += ipt_ssor_size[ i ] * param_int[nd][NEQ];
    else
      {
	E_Int       Ptomp = param_int[nd][PT_OMP];
	E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
	E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];
	E_Int cpt_actif   = 0;
	E_Int ithread_pre = 0;
	
	for(E_Int i = 0; i < ithread; i++)
	  if (param_int[nd][PtZoneomp + i] != - 2)
	    {
	      cpt_actif++;
	      if (cpt_actif > 1)
		ipt_ssor_shift += ipt_ssor_size[ ithread_pre ] * param_int[nd][NEQ];
	      ithread_pre = i;
	    }
      }
  }
