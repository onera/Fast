/*#include "Converter/IO/GenIO.h"
	   vector<E_Int> ni(nidom);  vector<E_Int> nj(nidom);   vector<E_Int> nk(nidom);
	   for (E_Int i = 0; i < nidom; i++) ni[i] = 1;
           vector<FldArrayF*> field(nidom);
	   vector<char*> zoneNames;
   	   for (E_Int i = 0; i < nidom; i++) field[i] = K_FLD::FldArrayF(size, nfld, iptro, true);
	   for (E_Int i = 0; i < nidom; i++) zonesNames = "zone";
           GenIO::tecwrite("toto.plt", "5e12", "ro,rou", field, zoneNmaes);
	   
*/
     

          ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
          nb_subzone    = ipt_nidom_loc [nitcfg-1];                                            //nbre sous-zone a la sousiter courante

            //if(ithread==1)printf("hum %d %d %d \n", nb_subzone, nd, nitcfg);
          //---------------------------------------------------------------------
          // -----Boucle sur param_int[nd][ ILES ] sous-zones ( //on skippe param_int[nd][ ILES ] parties qui converge + vite (Daude)
          // ---------------------------------------------------------------------
          for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
          {
            E_Int ndo   = nd;


            E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp;
            if (omp_mode == 1)
            { 
              E_Int       Ptomp = param_int[nd][PT_OMP];
              E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
              E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

              Nbre_thread_actif_loc = param_int[nd][ PtZoneomp  + Nbre_thread_actif ];
              ithread_loc           = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
              ipt_topo_omp          = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 1;
              ipt_inddm_omp         = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;


             //if (ithread == param_int[nd][IO_THREAD] && nitrun == 0 && ithread_loc != -1) printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
             //if (nd == 0 && nitrun == 0 && ithread_loc != -1) printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
             // printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
             //if (ithread == param_int[nd][IO_THREAD] && ithread_loc != -1) printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
              //ipt_inddm_omp[0],ipt_inddm_omp[1], ipt_inddm_omp[2],ipt_inddm_omp[3],ipt_inddm_omp[4],ipt_inddm_omp[5] );
             //if (ithread == 5 && ithread_loc != -1 && nitrun ==21)  printf("thread loc %d, Nthreads= %d, ithread= %d, zone= %d,  %d %d %d %d %d %d, iv= %d, jv= %d, nstep= %d, nitrun= %d \n",
             //ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
              //ipt_inddm_omp[0],ipt_inddm_omp[1], ipt_inddm_omp[2],ipt_inddm_omp[3],ipt_inddm_omp[4],ipt_inddm_omp[5], PtrIterOmp , PtZoneomp , nitcfg , PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6);

              if (ithread_loc == -1) {nd_current++; continue;}
            }

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

            //if(nd==2 && ithread==1) printf(" nd_ssz %d  %d \n",nd_subzone , nitcfg);
            //
            navier_stokes_struct_( ndo, nidom, Nbre_thread_actif_loc, ithread_loc, omp_mode, Nbre_socket, socket, mx_synchro , lssiter_verif, nptpsi, nitcfg, nitrun, first_it, nb_pulse, flagCellN,
                                  param_int[nd] , param_real[nd] ,
                                  temps               , ipt_tot       ,
                                  ipt_ijkv_sdm_thread , ipt_ind_dm_loc, ipt_ind_dm_socket, ipt_ind_dm_omp_thread, ipt_topology_socket, ipt_lok_thread, ipt_topo_omp, ipt_inddm_omp,
                                  ipt_cfl_thread      ,
                                  iptx[nd]                , ipty[nd]                , iptz[nd]            , iptCellN_loc     ,
                                  iptro[nd]               , iptro_m1[nd]            , iptrotmp[nd]        , iptro_ssiter[nd] ,
                                  iptmut[nd]              , 
                                  ipti[nd]                , iptj[nd]                , iptk[nd]            , iptvol[nd]       ,
                                  ipti0[nd]               , iptj0[nd]               , iptk0[nd]           , iptvol_df[nd]    ,
                                  iptventi[nd]            , iptventj[nd]            , iptventk[nd]        ,
                                  iptwig   + shift_wig    , iptstat_wig + shift_wig , iptrot+ shift_wig   ,
				  iptdrodm + shift_zone   , iptcoe  + shift_coe     , iptdelta[nd]        , iptro_res[nd]     );

            nd_current++;

            if(nd_current > mx_nidom)
             {
               if (ithread==1)
               {
                printf("------\n");
                printf("Error msg\n");
                printf("------\n");
                printf("resize MX_SSZONE. Present value= %d \n ", mx_nidom/nidom);
                printf("Value must be at least larger than : %d \n ", nd_current/nidom +2);
                printf("Just after the modules import of userscript.py, add the following python command:\n");
                printf("#\n");
                printf("#\n");
                printf("Fast.FastI.MX_SSZONE= %d\n ", nd_current/nidom +3);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }
             }
            //
          } //Fin boucle sur param_int[nd][ ILES ] sous-zones eventuelparam_int[nd][ ILES ] pour calcul RHS
