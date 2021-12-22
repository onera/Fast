/*#include "Converter/IO/GenIO.h"
	   vector<E_Int> ni(nidom);  vector<E_Int> nj(nidom);   vector<E_Int> nk(nidom);
	   for (E_Int i = 0; i < nidom; i++) ni[i] = 1;
           vector<FldArrayF*> field(nidom);
	   vector<char*> zoneNames;
   	   for (E_Int i = 0; i < nidom; i++) field[i] = K_FLD::FldArrayF(size, nfld, iptro, true);
	   for (E_Int i = 0; i < nidom; i++) zonesNames = "zone";
           GenIO::tecwrite("toto.plt", "5e12", "ro,rou", field, zoneNmaes);
	   
*/
     

          E_Int    cpu_perzone  =  nssiter*Nbre_thread_actif*2 + nd*(Nbre_thread_actif+1);
          E_Float* timer_omp_th = timer_omp + cpu_perzone + ithread;

          ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
          nb_subzone    = ipt_nidom_loc [nitcfg-1];                                            //nbre sous-zone a la sousiter courante
	 

           // cout <<  nb_subzone  << endl; 

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


              timer_omp_th = timer_omp + cpu_perzone + ithread_loc;

             //if (ithread == param_int[nd][IO_THREAD] && nitrun == 0 && ithread_loc != -1) printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
             //if (nd == 0 && nitrun == 0 && ithread_loc != -1) printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
             // printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
             //if (ithread == param_int[nd][IO_THREAD] && ithread_loc != -1) printf("thraed loxc %d %d %d %d %d %d %d %d %d %d \n", ithread_loc, Nbre_thread_actif_loc, ithread, nd, 
              //ipt_inddm_omp[0],ipt_inddm_omp[1], ipt_inddm_omp[2],ipt_inddm_omp[3],ipt_inddm_omp[4],ipt_inddm_omp[5] );
             //if (ithread == 24)  printf("thread loc %d, Nthreads= %d, ithread= %d, zone= %d %d ,  %d %d %d %d %d %d \n",
             //ithread_loc, Nbre_thread_actif_loc, ithread, nd, nd_subzone,
             // ipt_inddm_omp[0],ipt_inddm_omp[1], ipt_inddm_omp[2],ipt_inddm_omp[3],ipt_inddm_omp[4],ipt_inddm_omp[5]);

              if (ithread_loc == -1) { nd_current++; continue;}
            }


            //Init verrou rhs pour chaque sous zone et chaque thread actif:  init val to zero
            E_Int type = 4;
            E_Int* verrou_lhs_thread= verrou_lhs +             nd_current*Nbre_thread_actif + ithread_loc -1; 
            verrou_c_( verrou_lhs_thread, type);
            verrou_lhs_thread       = verrou_lhs + (mx_nidom + nd_current)*Nbre_thread_actif + ithread_loc -1; //pour calcul residu avant LU
            verrou_c_( verrou_lhs_thread, type );

            E_Int* ipt_lok_thread;
            //  Revoir cet adressage si scater et  socket>1 et ou nidom >1
            ipt_lok_thread   = ipt_lok   + nd_current*mx_synchro*Nbre_thread_actif;

            E_Int* ipt_ind_dm_loc         = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;      //ind_dm(6, < ssdom_lu,nssiter)
            E_Float* ipt_cfl_thread       = ipt_cfl         + (ithread_loc-1)*3+ ndo*3*Nbre_thread_actif;

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


	    //cout << "nitcfg= " << nitcfg << endl; 
	    //cout << ndo <<"  "<< ipt_ind_dm_omp_thread[0]<<" "<<ipt_ind_dm_omp_thread[1]<<" "<<ipt_ind_dm_omp_thread[2]<<" "<<ipt_ind_dm_omp_thread[3]<< endl;
	    //clock_t start, end;
	    
            //start = clock();
	    
            navier_stokes_struct_( ndo,    nidom, Nbre_thread_actif_loc, ithread_loc, ithread, omp_mode, layer_mode, Nbre_socket, socket, mx_synchro , 
                                   lssiter_verif, lexit_lu             ,nptpsi      , nitcfg , nssiter , nitrun    , first_it   , nb_pulse  , flagCellN,
                                  param_int[nd] , param_real[nd] ,
                                  temps               , ipt_tot       ,
                                  ipt_ijkv_sdm_thread , ipt_ind_dm_loc, ipt_ind_dm_socket, ipt_ind_dm_omp_thread, ipt_topology_socket, ipt_lok_thread, ipt_topo_omp, ipt_inddm_omp, timer_omp_th,
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


            //end = clock();
            //E_Float ns_duree;
//            E_Float ns_duree2;
//E_Float nb_pts = (ipt_ind_dm_omp_thread[1]+1)*(ipt_ind_dm_omp_thread[3]+1)*(ipt_ind_dm_omp_thread[5]+1);
//ns_duree  = (ns_end - ns_begin);///nb_pts;
            //cout <<"zone : "<< ndo <<" "<< "tps2 : "<< ns_duree2 << endl;

            //Flush Rhs
            E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX];
            //flush_real_( size , iptdrodm + shift_zone);
            if(nitcfg==1)
            {
              size = param_int[nd][NEQ_COE]*param_int[nd][NDIMDX];
              flush_real_( size , iptcoe + shift_coe);
            }
            //size = param_int[nd][NDIMDX];
            //flush_real_( size , iptmut[nd]);
            //#pragma omp flush
            //Go verrou rhs pour chaque sous zone et chaque thread actif: valeur mise a un
            type             = 1;
            verrou_lhs_thread= verrou_lhs + nd_current*Nbre_thread_actif + ithread_loc -1; 
            verrou_c_( verrou_lhs_thread, type );


            nd_current++;

            if(ithread_loc==1 && lexit_lu==0 && nitcfg*nitrun >15) {timer_omp[cpu_perzone]+=1; } //nbre echantillon

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
                printf("FastC.MX_SSZONE= %d\n ", nd_current/nidom +3);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }
             }
            //
          } //Fin boucle sur souszone pour calcul RHS
