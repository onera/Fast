          E_Int    cpu_perzone  =  nssiter*Nbre_thread_actif*2 + nd*(Nbre_thread_actif*2+1);
          E_Float* timer_omp_th = timer_omp + cpu_perzone + 1+ (ithread-1)*2;

          E_Int ndo   = nd;

          E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp;

            ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
            E_Int nd_subzone      = ipt_omp[ pttask + 1 ];
            Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
            ipt_topo_omp          = ipt_omp + pttask + 3 + Nbre_thread_actif ;
            ipt_inddm_omp         = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

            if (ithread_loc == -1) {continue;}
           
            ncells = (ipt_inddm_omp[1]-ipt_inddm_omp[0]+1)*(ipt_inddm_omp[3]-ipt_inddm_omp[2]+1)*(ipt_inddm_omp[5]-ipt_inddm_omp[4]+1);

            //if(nd==45 and nitrun >= 14401 ) printf("inddom %d %d %d %d %d %d %d %d ' ,nd_subz= ' %d %d %d %d %d %d %d %d \n", ipt_inddm_omp[0], ipt_inddm_omp[1],  ipt_inddm_omp[2], ipt_inddm_omp[3], ipt_inddm_omp[4], ipt_inddm_omp[5], ithread_loc, ithread, nd_subzone,Nbre_thread_actif_loc, nb_subzone, pttask, ptiter, ntask,nbtask,nitcfg );
            //printf("topo %d %d %d %d %d \n",ipt_topo_omp[0], ipt_topo_omp[1],  ipt_topo_omp[2], Nbre_thread_actif_loc, nd );
            //printf("shif %d %d %d  \n", shift_zone,shift_coe, shift_wig );

            //  Revoir cet adressage si scater et  socket>1 et ou nidom >1
            E_Int* ipt_lok_thread   = ipt_lok   + ntask*mx_synchro*Nbre_thread_actif;

            E_Int* ipt_ind_dm_loc         = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;      //ind_dm(6, < ssdom_lu,nssiter)
            E_Float* ipt_cfl_thread       = ipt_cfl         + (ithread_loc-1)*3+ ndo*3*Nbre_thread_actif;

            E_Float* iptCellN_loc; E_Int flagCellN;
            if (iptCellN[nd] == NULL) { flagCellN = 0; iptCellN_loc = iptro[nd];}
            else                      { flagCellN = 1; iptCellN_loc = iptCellN[nd]; }

            // Distribution de la sous-zone sur les threads
            indice_boucle_lu_(ndo, socket , Nbre_socket, lmin,
                              ipt_ind_dm_loc,
                              ipt_topology_socket, ipt_ind_dm_socket );

            navier_stokes_struct_( ndo,    Nbre_thread_actif_loc, ithread_loc, ithread, omp_mode, layer_mode, Nbre_socket, socket, mx_synchro , 
                                   lssiter_verif, lexit_lu             ,nptpsi      , nitcfg , nssiter , nitrun    , first_it   , nb_pulse  , flagCellN,
                                  param_int[nd] , param_real[nd] ,
                                  temps               ,
                                  ipt_ijkv_sdm_thread , ipt_ind_dm_loc, ipt_ind_dm_socket, ipt_inddm_omp, ipt_topology_socket, ipt_lok_thread, ipt_topo_omp, timer_omp_th,
                                  iptkrylov[nd]       , ipt_norm_kry[ithread-1],
                                  ipt_cfl_thread      ,
                                  iptx[nd]                , ipty[nd]                , iptz[nd]            , iptCellN_loc     , iptCellN_IBC[nd],
                                  iptro[nd]               , iptro_m1[nd]            , iptrotmp[nd]        , iptro_ssiter[nd] ,
                                  iptmut[nd]              , 
                                  ipti[nd]                , iptj[nd]                , iptk[nd]            , iptvol[nd]       ,
                                  ipti0[nd]               , iptj0[nd]               , iptk0[nd]           , iptvol_df[nd]    ,
                                  iptventi[nd]            , iptventj[nd]            , iptventk[nd]        ,
                                  iptwig   + shift_wig    , iptstat_wig + shift_wig , iptrot+ shift_wig   ,
				  iptdrodm + shift_zone   , iptcoe  + shift_coe     , iptdelta[nd]        , iptro_res[nd]  , iptsrc[nd]   );

            if(nitcfg==1)
            {
              E_Int size = param_int[nd][NEQ_COE]*param_int[nd][NDIMDX];
              flush_real_( size , iptcoe + shift_coe);
            }

            //Go verrou rhs pour chaque sous zone et chaque thread actif: valeur mise a un
            E_Int type             = 1;
            E_Int* verrou_lhs_thread= verrou_lhs + ntask*Nbre_thread_actif + ithread_loc -1; 
            verrou_c_( verrou_lhs_thread, type );

            if(ithread_loc==1 && lexit_lu==0 && nitcfg*nitrun >15 and (nitcfg < 3 or nitcfg == nssiter-1) ){ timer_omp[cpu_perzone]+=1; } //nbre echantillon

            if(ntask > mx_nidom)
             {
               if (ithread==1)
               {
                printf("------\n");
                printf("Error msg\n");
                printf("------\n");
                printf("resize MX_SSZONE. Present value= %d \n ", mx_nidom/nidom);
                printf("Value must be at least larger than : %d \n ", ntask/nidom +2);
                printf("Just after the modules import of userscript.py, add the following python command:\n");
                printf("#\n");
                printf("#\n");
                printf("FastC.MX_SSZONE= %d\n ", ntask/nidom +3);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }
             }
