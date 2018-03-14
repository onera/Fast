                 ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
                 nb_subzone     = ipt_nidom_loc [nitcfg-1];                                            //nbre sous-zone a la sousiter courante

                 E_Float* ipt_CL = iptro_CL[nd];
                 //---------------------------------------------------------------------
                 // -----Boucle sur param_int[nd][ ILES ] sous-zones ( //on skippe param_int[nd][ ILES ] parties qui converge + vite (Daude)
                 // ---------------------------------------------------------------------
                for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
                   {

                     E_Int* ipt_ind_dm_loc         = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;      //ind_dm(6, < ssdom_lu,nssiter)

                    // Distribution de la sous-zone sur les threads
                    indice_boucle_lu_(nd, ithread_loc, Nbre_thread_actif_loc, lmin,
                                      ipt_ind_dm_loc,
                                      ipt_topology_socket, ipt_ind_dm_omp_thread );

                    //if(ithread == param_int[nd][ IO_THREAD ]) { balance = balance +  (1+ipt_ind_dm_loc[1]-ipt_ind_dm_loc[0])*(1+ipt_ind_dm_loc[3]-ipt_ind_dm_loc[2])*(1+ipt_ind_dm_loc[5]-ipt_ind_dm_loc[4]);}

                    //if(nd==2 && ithread==1) printf(" dom %d %d %d %d %d %d \n", ipt_ind_dm_loc[0],ipt_ind_dm_loc[1],ipt_ind_dm_loc[2],ipt_ind_dm_loc[3],ipt_ind_dm_loc[4],ipt_ind_dm_loc[5]  );
                     //sortie de la carte d residu du Newton
                     //if( kimpli == 1)
                     // {

                         //if(ithread== param_int[nd][ IO_THREAD ]&& nd_subzone==nb_subzone-1  && nd==nidom-1) printf("balance %d %d \n", balance, nitcfg); 


                         if(lssiter_verif ==1  && nd_subzone ==0 && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
                          {
                            E_Int ijkv_lu[3];

                            ijkv_lu[0] = K_FUNC::E_max( 1, param_int[nd][ IJKV    ]/param_int[nd][ SIZE_SSDOM   ]);
                            ijkv_lu[1] = K_FUNC::E_max( 1, param_int[nd][ IJKV +1 ]/param_int[nd][ SIZE_SSDOM +1]);
                            ijkv_lu[2] = K_FUNC::E_max( 1, param_int[nd][ IJKV +2 ]/param_int[nd][ SIZE_SSDOM +2]);

                            E_Int ndim_rdm= ijkv_lu[0]*ijkv_lu[1]*ijkv_lu[2];
                            E_Int* ipt_nisdom_residu   =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter;                //nisdom_residu(nssiter)
                            E_Int* ipt_it_bloc         =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter*2;    //it_bloc(nidom)

                            E_Int* ipt_it_lu_ssdom_loc =  ipt_it_lu_ssdom[nd];
                            E_Int* ipt_it_target_ssdom =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ];
                            E_Int* ipt_it_target_old   =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ]*2;
                            E_Int* ipt_it_temp_ssdom   =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ]*3;
                            E_Int* ipt_no_lu           =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ]*4;

			    //                           if(ithread== param_int[nd][ IO_THREAD ] && nd_subzone==nb_subzone-1  && nd==nidom-1) printf("balance %d %d \n", balance, nitcfg); 
 
                            cprdu3s1_(nd,nitcfg, nssiter          , param_int[nd][ NEQ ], param_int[nd][ NDIMDX ], ndim_rdm , nitrun,
                                      param_int[nd][ IFLOW ]      , param_int[nd][ ILES ] , ithread_loc          , omp_mode,
                                      param_int[nd][ MXSSDOM_LU ] , ipt_it_bloc[0]        , param_real[nd][ EPSI_NEWTON ] ,
                                      param_int[nd]+ SIZE_SSDOM   , ipt_nisdom_residu     ,
                                      param_int[nd]+ NIJK         , param_int[nd]+ IJKV     ,ijkv_lu     ,
                                      ipt_ind_dm_loc       ,
                                      ipt_it_lu_ssdom_loc  , ipt_it_target_ssdom  , ipt_it_target_old     ,
                                      ipt_it_temp_ssdom    , ipt_no_lu            ,
                                      iptdrodm + shift_zone   , iptrdm[nd]); 
                          }

                     if( kimpli == 1 )
                      {
                         //
                         //
                         // CL sur rhs pour implicitation
                         E_Int lrhs=1; E_Int lcorner=1; E_Int ipt_ind_coe_thread[6];
                         BCzone( nd, lrhs, lcorner,
                                 param_int[nd], param_real[nd],
                                 npass,
                                 ipt_ind_dm_loc         , ipt_ind_dm_omp_thread      ,
                                 ipt_ind_CL_thread      , ipt_ind_CL119_thread       , ipt_ind_coe_thread,
                                 iptdrodm + shift_zone  , ipti[nd]                   , iptj[nd]                 , iptk[nd]       ,
                                 iptx[nd]               , ipty[nd]                   , iptz[nd]                 ,
                                 iptventi[nd]           , iptventj[nd]               , iptventk[nd]             );

                         if(lcorner  == 0 )correct_coins_(nd, param_int[nd], ipt_ind_coe_thread , iptdrodm + shift_zone );

                         if(lexit_lu == 0 ) invlu_(nd                     , nitcfg      ,nitrun, param_int[nd], param_real[nd],
                                                   ipt_ind_coe_thread     ,
                                                   iptrotmp[nd]           , iptro_ssiter[nd]        , iptdrodm + shift_zone ,
                                                   ipti[nd]               , iptj[nd]                , iptk[nd]              ,
                                                   iptventi[nd]           , iptventj[nd]            , iptventk[nd]          ,
                                                   iptcoe  + shift_coe   );
                      } //fin kimpli


                     //
                     //
                     // CL sur var primitive
/*                     E_Int lrhs=0; E_Int lcorner=0;
                     BCzone( nd, lrhs, lcorner,
                             param_int[nd], param_real[nd],
                             npass,
                             ipt_ind_dm_loc         , ipt_ind_dm_omp_thread      ,
                             ipt_ind_CL_thread      , ipt_ind_CL119_thread       ,  ipt_ind_coe_thread,
                             ipt_CL                 , ipti[nd]                   , iptj[nd]                 , iptk[nd]       ,
                             iptx[nd]               , ipty[nd]                   , iptz[nd]                 ,
                             iptventi[nd]           , iptventj[nd]               , iptventk[nd]             );
*/
                     // Selective Frequency Damping
                     if(( (nitcfg == nssiter && lssiter_verif==1) || (nitcfg == nssiter-1 && lssiter_verif==0)) && param_int[nd][SFD] == 1)
                     {
                       sfd_(param_int[nd], param_real[nd], nitrun, ipt_ind_dm_omp_thread, ipt_CL, iptrof[nd], iptcoe + shift_coe, iptvol[nd]);
                     }

                   } //fin boucle sous-zone
