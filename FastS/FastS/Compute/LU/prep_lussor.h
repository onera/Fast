if (param_int[ nd ][ NB_RELAX ] > 1)
  {
    ipt_ssor_shift = iptssor[ nd ];
    ipt_ssortmp_shift = iptssortmp[ nd ];
    E_Int indice; E_Int shift=0;

    if (omp_mode == 0)
      for (E_Int i = 0; i < ithread - 1; i++)
	{
	  indice = nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
	  ipt_ssor_shift += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
	  ipt_ssortmp_shift += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
	  shift += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
	}
    else
      {
	E_Int       Ptomp = param_int[ nd ][ PT_OMP ];
	E_Int  PtrIterOmp = param_int[ nd ][ Ptomp + nitcfg - 1 ];   
	E_Int  PtZoneomp  = param_int[ nd ][ PtrIterOmp + nd_subzone ];
	E_Int cpt_actif   = 0;
	E_Int ithread_pre = 0;
	
	for(E_Int i = 0; i < ithread; i++)
	  if (param_int[ nd ][ PtZoneomp + i ] != - 2)
	    {
	      cpt_actif++;
	      if (cpt_actif > 1)
		{
		  indice = nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + ithread_pre;
		  ipt_ssor_shift += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
		  ipt_ssortmp_shift += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
		}
	      ithread_pre = i;
	    }
      }

    if (param_int[ nd ][ LU_MATCH ] == 1)
      {
#pragma omp barrier
      }
    
    ssor_size = ipt_ssor_size[ nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + ithread - 1 ];

    //printf("shift = %d et thread = %d \n", shift/5,ithread);
  }

 else
   ssor_size = 1;
