      //sortie de la carte d residu du Newton
      if(lssiter_verif ==1  && nd_subzone ==0 && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
       {
        if (nb_subzone !=1)// pour 1 drodm OK pour calcul residu
        {
          E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX];
          flush_real_( size , iptdrodm + shift_zone);
        }

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


         //protection pour effet de bord OMP quand grand nombre de thread test souvent le meme verrou
         if(ndim_rdm < Nbre_thread_actif_loc){barrier_residu=1;}
         
         if (ithread_loc != -1)
         { 
           //rhs_begin = omp_get_wtime();

           cprdu3s1_(nd,nitcfg, nssiter          , param_int[nd][ NEQ ], param_int[nd][ NDIMDX ], ndim_rdm , nitrun,
                   param_int[nd][ IFLOW ]      , param_int[nd][ ILES ] ,
                   ithread_loc                 , Nbre_thread_actif_loc , omp_mode,
                   param_int[nd][ MXSSDOM_LU ] , ipt_it_bloc[0]        , param_real[nd][ EPSI_NEWTON ] ,
                   param_int[nd]+ SIZE_SSDOM   , ipt_nisdom_residu     ,
                   param_int[nd]+ NIJK         , param_int[nd]+ IJKV   , ijkv_lu     ,
                   ipt_ind_dm_loc       ,
                   ipt_it_lu_ssdom_loc  , ipt_it_target_ssdom  , ipt_it_target_old     ,
                   ipt_it_temp_ssdom    , ipt_no_lu            ,
                   iptdrodm + shift_zone   , iptrdm[nd]); 

           //rhs_end = omp_get_wtime();
           //if(ithread_loc == 1) {printf(" time residu= %g  %d  %d \n",(rhs_end - rhs_begin), nitcfg, ithread);}

           //Go verrou residu pour chaque sous zone et chaque thread actif pour ne pas attaquer LU avant fin calcul residu en mode1: 
           E_Int type   = 1;
           E_Int* verrou_lhs_thread= verrou_lhs +(nbtask + ntask)*Nbre_thread_actif + ithread_loc -1; verrou_c_( verrou_lhs_thread, type );

         }// test ithread_loc =-1
       }
