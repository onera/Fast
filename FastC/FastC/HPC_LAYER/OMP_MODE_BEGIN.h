        nbtask = ipt_omp[nitcfg-1]; 
        ptiter = ipt_omp[nssiter+ nitcfg-1];

        for (E_Int ntask = 0; ntask < nbtask; ntask++)
          {
             E_Int pttask     = ptiter + ntask*(6+Nbre_thread_actif*7);
             E_Int nd         = ipt_omp[ pttask ];
             E_Int nd_subzone = ipt_omp[ pttask + 1 ];

             E_Int* ipt_ind_dm_loc = ipt_ind_dm[nd] + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
             E_Int* ipt_nidom_loc  = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
             E_Int  nb_subzone     = ipt_nidom_loc [nitcfg-1]; 

             E_Float* ipt_ssor_shift; E_Float* ipt_ssortmp_shift; E_Int ssor_size; 
             E_Int* ipt_ind_dm_thread;

             E_Int lmin = 10;
             if (param_int[nd][ITYPCP] == 2) lmin = 4;

             shift_zone=0; shift_coe=0; shift_mu=0; shift_wig=0;
             for (E_Int n = 0; n < nd; n++)
             {
              shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
              shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
              shift_mu   = shift_mu   + param_int[n][ NDIMDX ];
              if(param_int[n][ KFLUDOM ]==2){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*3;}
             }


             ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
             Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
             ipt_ind_dm_thread     = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

             if (ithread_loc == -1) {continue;}

