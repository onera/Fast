 E_Int indice;
if (param_int[ nd ][ NB_RELAX ] > 1 || param_int[ nd ][ LU_MATCH ]==1 )
  {
    ipt_ssor_shift    = iptssor[    nd ];
    ipt_ssortmp_shift = iptssortmp[ nd ];

    for (E_Int i = 0; i < ithread_loc - 1; i++)
     {
	  indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
	  ipt_ssor_shift    += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
	  ipt_ssortmp_shift += ipt_ssor_size[ indice ] * param_int[ nd ][ NEQ ];
     }

    if (param_int[ nd ][ LU_MATCH ] == 1)
      {
#pragma omp barrier
      }
    
    indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + ithread_loc -1;

    ssor_size = ipt_ssor_size[ indice ];

    //printf("shift = %d et thread = %d \n", shift/5,ithread);
  }

 else
   ssor_size = param_int[nd][NDIMDX];
