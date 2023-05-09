                ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
                Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
                ipt_topo_omp          = ipt_omp + pttask + 3 + Nbre_thread_actif ;
                ipt_inddm_omp         = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

                if (ithread_loc == -1) {continue;}

                E_Int* ind_mtr  = ipt_ind_sdm     + (ithread-1)*6;
                for (E_Int i = 0; i < 6; i++) ind_mtr[i] = ipt_inddm_omp[i];

                if (ind_mtr[0]==1) ind_mtr[0] -= param_int[nd][NIJK+3];
                if (ind_mtr[2]==1) ind_mtr[2] -= param_int[nd][NIJK+3];
                if (ind_mtr[4]==1) ind_mtr[4] -= param_int[nd][NIJK+4];
                if (ind_mtr[1]== param_int[nd][IJKV  ]) ind_mtr[1] += param_int[nd][NIJK+3];
                if (ind_mtr[3]== param_int[nd][IJKV+1]) ind_mtr[3] += param_int[nd][NIJK+3];
                if (ind_mtr[5]== param_int[nd][IJKV+2]) ind_mtr[5] += param_int[nd][NIJK+4];

                //tableau metric d'epaisseur k =1 si 3Dhomogene, cartesien ou 2D
                if      (param_int[nd][ITYPZONE ] == 1 )
                {    ind_mtr[4]= 1; ind_mtr[5]= 1;
                     E_Int ssblk_ktranche = ipt_topo_omp[0]*ipt_topo_omp[1];//seul les threads du bloc K=1 calculent
                     if (ithread_loc > ssblk_ktranche ) continue; 
                }
                else if(param_int[nd][ ITYPZONE ] == 2 )
                { 
                   for (E_Int i = 0; i < 6; i++) ind_mtr[i] =1;
                   if ( ithread_loc !=1 ) continue; 
                }
