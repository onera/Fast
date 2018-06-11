                      E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
                      E_Int  nb_subzone    = ipt_nidom_loc [nitcfg-1];                                           //nbre sous-zone a la sousiter courante
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
                         }// omp_mode
