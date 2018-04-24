                E_Int shift_omp = ipt_param_int[nd][ PT_OMP ];

                Nbre_thread_actif_loc = ipt_param_int[nd][ shift_omp  + Nbre_thread_actif ];
                ithread_loc           = ipt_param_int[nd][ shift_omp  +  ithread -1       ] +1 ;
                ipt_topo_omp          = ipt_param_int[nd] + shift_omp +  Nbre_thread_actif + 1;
                ipt_inddm_omp         = ipt_param_int[nd] + shift_omp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;
                if (ithread_loc == -1) {continue;}

                E_Int* ind_mtr  = ipt_ind_sdm     + (ithread-1)*6;
                for (E_Int i = 0; i < 6; i++) ind_mtr[i] = ipt_inddm_omp[i];

                if (ind_mtr[0]==1) ind_mtr[0] -= ipt_param_int[nd][NIJK+3];
                if (ind_mtr[2]==1) ind_mtr[2] -= ipt_param_int[nd][NIJK+3];
                if (ind_mtr[4]==1) ind_mtr[4] -= ipt_param_int[nd][NIJK+4];
                if (ind_mtr[1]== ipt_param_int[nd][IJKV  ]) ind_mtr[1] += ipt_param_int[nd][NIJK+3];
                if (ind_mtr[3]== ipt_param_int[nd][IJKV+1]) ind_mtr[3] += ipt_param_int[nd][NIJK+3];
                if (ind_mtr[5]== ipt_param_int[nd][IJKV+2]) ind_mtr[5] += ipt_param_int[nd][NIJK+4];

                //tableau metric d'epaisseur k =1 si 3Dhomogene, cartesien ou 2D
                if      (ipt_param_int[nd][ITYPZONE ] == 1 )
                {    ind_mtr[4]= 1; ind_mtr[5]= 1;
                     E_Int ssblk_ktranche = ipt_topo_omp[0]*ipt_topo_omp[1];//seul les thraed du bloc K=1 calculent
                     if (ithread_loc > ssblk_ktranche ) continue; 
                }
                else if(ipt_param_int[nd][ ITYPZONE ] == 2 )
                { 
                   for (E_Int i = 0; i < 6; i++) ind_mtr[i] =1;
                   if ( ithread_loc !=1 ) continue; 
                }
