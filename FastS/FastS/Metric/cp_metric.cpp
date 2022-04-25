        E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp; 
        E_Int* ipt_ind_sdm  = ipt_ind_CL;

        E_Int ipt_ind_dm_loc[6];           
        E_Int ipt_topology_socket_thread[6];
        E_Int ipt_ind_sdm_thread[6];
        E_Int ipt_ind_coe_thread[6];
        E_Int ipt_ind_grad_thread[6];

        if (omp_mode == 1)
            { 
              // loop calcul normale
              E_Int barrier = 0;
              for (E_Int nd = 0; nd < nidom; nd++)
              {
               if(param_int[nd][LALE]==2 && param_int[nd][ITYPZONE]!=4)
               {
                 barrier = 1;
#                include "FastS/Metric/indice_omp1.h" 
                 cp_tijk_( param_int[nd], iptx[nd], ipty[nd], iptz[nd], ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], ind_mtr);

               }
               if (nd == nidom-1 && barrier == 1)
               {
       	        #pragma omp barrier
               }
              }
              if( barrier == 1)
              {
       	       #pragma omp barrier
              }
              // loop calcul volume
              barrier = 0;
              for (E_Int nd = 0; nd < nidom; nd++)
              {
               if(param_int[nd][LALE]==2 && param_int[nd][ITYPZONE]!=4)
               {
                 barrier = 1;
#                include "FastS/Metric/indice_omp1.h" 
                 cp_vol_(  param_int[nd], iptx[nd], ipty[nd], iptz[nd], ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], iptvol[nd], ind_mtr);

                 for (E_Int k = ind_mtr[4]; k <= ind_mtr[5]; k++){ 
                  for (E_Int j = ind_mtr[2]; j <= ind_mtr[3]; j++){ 
                   for (E_Int i = ind_mtr[0]; i <= ind_mtr[1]; i++){ 

                     E_Int l =  (i+ param_int[nd][NIJK_MTR+3]-1)*param_int[nd][NIJK_MTR]
                              + (j+ param_int[nd][NIJK_MTR+3]-1)*param_int[nd][NIJK_MTR+1]
                              + (k+ param_int[nd][NIJK_MTR+4]-1)*param_int[nd][NIJK_MTR+2];

                     iptvol[nd][l] = K_FUNC::E_max(iptvol[nd][l], 1.e-30);
                    }
                   }
                  }
               }//ale
               if(nd == nidom-1 && barrier == 1)
               {
       	         #pragma omp barrier
               }
              }//loop zone
              if(barrier ==1)
              {
       	         #pragma omp barrier
              }

              // loop extrapolation
              for (E_Int nd = 0; nd < nidom; nd++)
              {
               if(param_int[nd][LALE]==2 && param_int[nd][ITYPZONE]!=4)
               {
#               include "FastS/Metric/indice_omp1.h" 
                if (param_int[nd][ ITYPZONE ] != 2  && ithread_loc == 1)
                {
                 ipt_ind_dm_loc[0]   = 1; 
                 ipt_ind_dm_loc[2]   = 1; 
                 ipt_ind_dm_loc[4]   = 1; 
                 if(param_int[nd][ ITYPZONE ] == 0) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = param_int[nd][ IJKV+2 ];
                   }
                 else if(param_int[nd][ ITYPZONE ] == 1) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1;
                   }
                 else if(param_int[nd][ ITYPZONE ] == 2) 
                   {
                    ipt_ind_dm_loc[1]   = 1 ;
                    ipt_ind_dm_loc[3]   = 1; 
                    ipt_ind_dm_loc[5]   = 1; 
                   }
                 else
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1 ;
                   }

                 E_Int* ipt_nijk_xyz = param_int[nd]+ NIJK_XYZ;
                 E_Int* ipt_nijk_mtr = param_int[nd]+ NIJK_MTR;
                 E_Int* ipt_nijk     = param_int[nd]+ NIJK;

                 tijk_extrap_( param_int[nd][ NDIMDX_MTR ], param_int[nd][ NDIMDX_XYZ ] , ipt_nijk_xyz, ipt_nijk_mtr,
                               param_int[nd][ NEQ_IJ ]    , param_int[nd][ NEQ_K ],
                               ipt_ind_dm_loc,
                               ipt_degen[nd] ,
                               ipti[nd]      , iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], iptvol[nd]); 

                 if(param_int[nd][ IFLOW ] ==3)
                   { ipt_ind_dm_loc[1]= param_int[nd][ IJKV ]; ipt_ind_dm_loc[3]= param_int[nd][ IJKV+1 ]; ipt_ind_dm_loc[5]= param_int[nd][IJKV+2];

                     dist_extrap_( param_int[nd][ NDIMDX ], param_int[nd][ NDIMDX_XYZ ] , ipt_nijk, ipt_nijk_xyz,
                                    ipt_ind_dm_loc, ipt_degen[nd] , iptmut[nd]+param_int[nd][NDIMDX]);
                   }
                }
               }//ale
              }//loop zone
            }
        else
            { //omp_mode 0
              for (E_Int nd = 0; nd < nidom; nd++)
              {
               if(param_int[nd][LALE]==2 && param_int[nd][ITYPZONE]!=4)
               {
                E_Int lmin = 10;
                if (param_int[nd][ ITYPCP ] == 2) lmin = 4;

                ipt_ind_dm_loc[0]   = 1; 
                ipt_ind_dm_loc[2]   = 1; 
                ipt_ind_dm_loc[4]   = 1; 
                if(param_int[nd][ ITYPZONE ] == 0) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = param_int[nd][ IJKV+2 ];
                   }
                else if(param_int[nd][ ITYPZONE ] == 1) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1;
                   }
                else if(param_int[nd][ ITYPZONE ] == 2) 
                   {
                    ipt_ind_dm_loc[1]   = 1 ;
                    ipt_ind_dm_loc[3]   = 1; 
                    ipt_ind_dm_loc[5]   = 1; 
                   }
                else
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1 ;
                   }

                 indice_boucle_lu_(nd, socket , Nbre_socket, lmin,
                                   ipt_ind_dm_loc, 
                                   ipt_topology_socket_thread, ipt_ind_dm_socket );

                 skmtr_( nd, param_int[nd], param_real[nd],
	                iptx[nd], ipty[nd], iptz[nd], ipt_degen[nd], iptmut[nd]+param_int[nd][NDIMDX],
                        ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], iptvol[nd], iptventi[nd], iptventj[nd], iptventk[nd],
                        ipt_ijkv_sdm_thread,
                        ipt_ind_sdm_thread , ipt_ind_coe_thread, ipt_ind_grad_thread    , 
                        ipt_ind_dm_loc     , ipt_ind_dm_socket , ipt_ind_dm_omp_thread  ,
                        ipt_topology_socket_thread, 
                        ithread_sock       , thread_parsock    , Nbre_thread_actif, Nbre_socket, socket,
                        ithread);

                 // cutoff vol
                 #pragma omp for
                 for (E_Int i = 0; i < param_int[nd][NDIMDX_MTR]; i++) { iptvol[nd][i] = K_FUNC::E_max(iptvol[nd][i], 1.e-30);}

                }//ale
              }//zone
           }// omp_mode


