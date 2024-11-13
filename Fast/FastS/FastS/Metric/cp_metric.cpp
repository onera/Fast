        E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp; 
        E_Int* ipt_ind_sdm   = ipt_ind_CL;

        // loop calcul normale
        for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
          E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
          E_Int nd = ipt_omp[ pttask ];

           if(param_int[nd][LALE]>=2 && param_int[nd][ITYPZONE]!=4)
           {
#             include "FastC/Metric/indice_omp1.h" 
             cp_tijk_( param_int[nd], iptx[nd], ipty[nd], iptz[nd], ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], ind_mtr);
           }
        }
        E_Int barrier_tijk = 0; E_Int barrier_vol = 0;
        for (E_Int nd = 0; nd < nidom; nd++)
          { if(param_int[nd][LALE]>=2 && param_int[nd][ITYPZONE]!=4) { barrier_tijk = 1;}
            if(param_int[nd][LALE]==3 && param_int[nd][ITYPZONE]!=4) { barrier_vol  = 1;}
          }

        if( barrier_tijk == 1)
        {
       	 #pragma omp barrier
        }
        // loop calcul volume
        for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
          E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
          E_Int nd = ipt_omp[ pttask ];

           if(param_int[nd][LALE]==3 && param_int[nd][ITYPZONE]!=4)
           {
             //recuperation pointeur pour stockage volume instant N+1
             E_Float* vol_p = iptvol[nd] + param_int[nd][PT_VOL]* param_int[nd][ NDIMDX_MTR ];

#            include "FastC/Metric/indice_omp1.h" 
             cp_vol_(  param_int[nd], iptx[nd], ipty[nd], iptz[nd], ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], vol_p, ind_mtr);

           }//ale
        }//loop task 
        if(barrier_vol ==1)
        {
       	   #pragma omp barrier
        }

        // loop extrapolation
        // loop calcul volume
        for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
          E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
          E_Int nd = ipt_omp[ pttask ];
              
          E_Int* ipt_ind_dm_loc= shift_lu + 6*nd*Nbre_thread_actif + (ithread-1)*6;

         if(param_int[nd][LALE]>=2 && param_int[nd][ITYPZONE]!=4)
         {
          E_Float* vol_p = iptvol[nd] + param_int[nd][PT_VOL]* param_int[nd][ NDIMDX_MTR ];
#         include "FastC/Metric/indice_omp1.h" 
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

           tijk_extrap_( param_int[nd][ NDIMDX_MTR ], param_int[nd][ NDIMDX_XYZ ] , ipt_nijk_xyz, ipt_nijk_mtr, ipt_nijk,
                         param_int[nd][ NEQ_IJ ]    , param_int[nd][ NEQ_K ],
                         ipt_ind_dm_loc,
                         ipt_degen[nd] ,
                         ipti[nd]      , iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], vol_p); 

           if(param_int[nd][ IFLOW ] ==3 && param_int[nd][LALE]==3)
             { ipt_ind_dm_loc[1]= param_int[nd][ IJKV ]; ipt_ind_dm_loc[3]= param_int[nd][ IJKV+1 ]; ipt_ind_dm_loc[5]= param_int[nd][IJKV+2];

               dist_extrap_( param_int[nd][ NDIMDX ], param_int[nd][ NDIMDX_XYZ ] , ipt_nijk, ipt_nijk_xyz,
                              ipt_ind_dm_loc, ipt_degen[nd] , iptmut[nd]+param_int[nd][NDIMDX]);
             }
          }
         }//ale
        }//loop task 
        if(barrier_tijk ==1)
        {
       	   #pragma omp barrier
        }

