          ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
          nb_subzone    = ipt_nidom_loc [nitcfg-1];                                            //nbre sous-zone a la sousiter courante

          //---------------------------------------------------------------------
          // -----Boucle sur param_int[nd][ ILES ] sous-zones ( //on skippe param_int[nd][ ILES ] parties qui converge + vite (Daude)
          // ---------------------------------------------------------------------
          for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
          {
            E_Int ndo   = nd;


            E_Int* ipt_lok_thread;
            //  Revoir cet adressage si scater et  socket>1 et ou nidom >1
            ipt_lok_thread   = ipt_lok   + nd_current*mx_synchro*Nbre_thread_actif;

            E_Int* ipt_ind_dm_loc         = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;      //ind_dm(6, < ssdom_lu,nssiter)
            E_Float* ipt_cfl_thread       = ipt_cfl         + (ithread-1)*3+ ndo*3*Nbre_thread_actif;

            E_Float* iptCellN_loc; E_Int flagCellN;
            if (iptCellN[nd] == NULL) { flagCellN = 0; iptCellN_loc = iptro[nd];}
            else                      { flagCellN = 1; iptCellN_loc = iptCellN[nd]; }

            // Distribution de la sous-zone sur les threads
            //E_Int icp_loc =2;
            
            //Pour les CL
            indice_boucle_lu_(ndo, ithread_loc, Nbre_thread_actif_loc, lmin,
                              ipt_ind_dm_loc,
                              ipt_topology_socket, ipt_ind_dm_omp_thread );

            indice_boucle_lu_(ndo, socket , Nbre_socket, lmin,
                              ipt_ind_dm_loc,
                              ipt_topology_socket, ipt_ind_dm_socket );

                     // CL sur var primitive
/*                     E_Int lrhs=0; E_Int lcorner=0;E_Int npass_loc =0;
                     E_Float* ipt_CL = iptro_CL[nd];


                     BCzone( nd, lrhs, lcorner,
                             param_int[nd], param_real[nd],
                             npass_loc,
                             ipt_ind_dm_loc         , ipt_ind_dm_omp_thread      ,
                             ipt_ind_CL_thread      , ipt_ind_CL119_thread       , ipt_ind_coe_thread,
                             iptro_ssiter[nd]       , ipti[nd]                   , iptj[nd]                 , iptk[nd]       ,
                             //iptro_CL[nd]       , ipti[nd]                   , iptj[nd]                 , iptk[nd]       ,
                             iptx[nd]               , ipty[nd]                   , iptz[nd]                 ,
                             iptventi[nd]           , iptventj[nd]               , iptventk[nd]             );
*/
            //if(nd==2 && ithread==1) printf(" nd_ssz %d  %d \n",nd_subzone , nitcfg);

            navier_stokes_struct_( ndo, nidom, Nbre_thread_actif_loc, ithread_loc, Nbre_socket, socket, mx_synchro , lssiter_verif, nptpsi, nitcfg, nitrun, first_it, nb_pulse, flagCellN,
                                  param_int[nd] , param_real[nd] ,
                                  temps               , ipt_tot                 ,
                                  ipt_ijkv_sdm_thread , ipt_ind_dm_loc      , ipt_ind_dm_socket       , ipt_ind_dm_omp_thread  ,  ipt_topology_socket , ipt_lok_thread       ,
                                  ipt_cfl_thread      ,
                                  iptx[nd]                , ipty[nd]                , iptz[nd]            , iptCellN_loc     ,
                                  iptro[nd]               , iptro_m1[nd]            , iptrotmp[nd]        , iptro_ssiter[nd] ,
                                  iptmut[nd]              , iptdist[nd]             ,
                                  ipti[nd]                , iptj[nd]                , iptk[nd]            , iptvol[nd]       ,
                                  ipti0[nd]               , iptj0[nd]               , iptk0[nd]           , iptvol_df[nd]    ,
                                  iptventi[nd]            , iptventj[nd]            , iptventk[nd]        ,
                                  iptwig   + shift_wig    , iptstat_wig + shift_wig , iptrot+ shift_wig   ,
				  iptdrodm + shift_zone   , iptcoe  + shift_coe     );

            nd_current++;

            if(nd_current > mx_nidom) {printf("redimensionner mx_nidom= %d a %d\n ", mx_nidom ,nd_current); exit(0);}
            //
          } //Fin boucle sur param_int[nd][ ILES ] sous-zones eventuelparam_int[nd][ ILES ] pour calcul RHS
