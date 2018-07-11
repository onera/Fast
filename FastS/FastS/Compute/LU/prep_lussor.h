if (param_int[ nd ][ NB_RELAX ] > 1)
  {
    ipt_ssor_shift = iptssor[nd];

    if (omp_mode == 0)
      for (E_Int i = 0; i < ithread - 1; i++)
	{
	  indice = nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
	  ipt_ssor_shift += ipt_ssor_size[ indice ] * param_int[nd][NEQ];
	}
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
		{
		  indice = nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + ithread_pre;
		  ipt_ssor_shift += ipt_ssor_size[ indice ] * param_int[nd][NEQ];
		}
	      ithread_pre = i;
	    }
      }
  }
