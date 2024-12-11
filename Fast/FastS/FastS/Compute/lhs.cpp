if( kimpli == 1  && param_int[0][LU_MATCH]==1 && param_int_tc != NULL)
  { 
   #pragma omp master
    { //Raccord V0
      
      E_Int rk           = param_int[0][RK];
      E_Int exploc       = param_int[0][EXPLOC];
      E_Int numpassage   = 1;

      K_FASTC::setInterpTransfersFast(iptdrodm_transfer, vartype, param_int_tc,
                         param_real_tc, param_int, param_real, ipt_omp, linelets_int, linelets_real,
      			 it_target, nidom, ipt_timecount, mpi, nitcfg , nssiter, rk, exploc, numpassage);
    }
    #pragma omp barrier
  }

    E_Float* ipt_ssor_shift; E_Float* ipt_ssortmp_shift; E_Int ssor_size;
    //calcul residu si necessaire
    E_Int barrier_residu = 0; 
       for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
         E_Int pttask     = ptiter + ntask*(6+Nbre_thread_actif*7);
         E_Int nd         = ipt_omp[ pttask ];
         E_Int nd_subzone = ipt_omp[ pttask + 1 ];

         E_Int* ipt_ind_dm_loc = ipt_ind_dm[nd] + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
         E_Int* ipt_nidom_loc  = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
         E_Int  nb_subzone     = ipt_nidom_loc [nitcfg-1]; 

         shift_zone=0; shift_coe=0;
         for (E_Int n = 0; n < nd; n++)
          {
           shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
           shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
          }

         if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
          {
#ifdef _OPENMP
            E_Float lhs_begin = omp_get_wtime();
#else
            E_Float lhs_begin = 0.;
#endif
            ithread_loc              = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
            Nbre_thread_actif_loc    = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
            E_Int* ipt_ind_dm_thread = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

            if(nitcfg==1 && ithread_loc != -1) { E_Int size =param_int[nd][NEQ_COE]*param_int[nd][NDIMDX]; flush_real_( size, iptcoe + shift_coe); }

            if(ithread_loc != -1)
             {
              // CL sur rhs pour implicitation
              if( kimpli == 1 )
               {  
                 E_Int lrhs=1; E_Int lcorner=1;
                 E_Int* ipt_shift_lu   = shift_lu + 6*ntask*Nbre_thread_actif + (ithread-1)*6;
                 K_FASTS::BCzone( nd, lrhs, nitcfg, lcorner, 
                                  param_int[nd], param_real[nd],
                                  npass,
                                  ipt_ind_dm_loc         , ipt_ind_dm_thread      ,
                                  ipt_ind_CL_thread      , ipt_ind_CL119          , ipt_ind_CLgmres, ipt_shift_lu,
                                  iptdrodm + shift_zone  , ipti[nd]               , iptj[nd]       , iptk[nd]          ,
                                  iptx[nd]               , ipty[nd]               , iptz[nd]       ,
                                  iptventi[nd]           , iptventj[nd]           , iptventk[nd]   , iptrotmp[nd], iptmut[nd]);  

                 //for (E_Int n = 0; n < 6 ;n++){ ipt_shift_lu[n]= ipt_ind_dm_thread[n];}
                 if(lcorner  == 0 )correct_coins_(nd, param_int[nd], ipt_shift_lu , iptdrodm + shift_zone );
                 //if(lcorner  == 0 )correct_coins_(nd, param_int[nd], ipt_ind_dm_thread , iptdrodm + shift_zone );

               }//kimpli
             }//ithread_loc
          }// test skip lbm et unstruct

         if(lssiter_verif ==1 )
          {
            if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4) // on skippe les eventuelles zone non structurees ou LBM
            {
#include   "FastS/Compute/verrou_residus.h"
//            if( ntask==0 )
//             {
//#pragma omp barrier
//             }
               E_Int* ipt_topo_omp; E_Int* ipt_ind_dm_thread;

               ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
               Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
               ipt_ind_dm_thread     = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

              //sortie de la carte residu du Newton
#include      "FastS/Compute/residus_navier.h"
            } // test lbm/unstructured
          }// test residu

        }// loop task residu

    if(barrier_residu==1)
    {
     #pragma omp barrier
    }
    
    for (E_Int ntask = 0; ntask < nbtask; ntask++)
      {
        
        E_Int pttask     = ptiter + ntask*(6+Nbre_thread_actif*7);
        E_Int nd         = ipt_omp[ pttask ];
        E_Int nd_subzone = ipt_omp[ pttask + 1 ];

        E_Int* ipt_ind_dm_loc = ipt_ind_dm[nd] + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
        E_Int* ipt_nidom_loc  = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
        E_Int  nb_subzone     = ipt_nidom_loc [nitcfg-1]; 

        E_Int lmin = 10;
        if (param_int[nd][ITYPCP] == 2) lmin = 4;

        shift_zone=0; shift_coe=0;
        for (E_Int n = 0; n < nd; n++)
        {
         shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
         shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
        }

        if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
        {
#ifdef _OPENMP
            E_Float lhs_begin = omp_get_wtime();
#else
            E_Float lhs_begin = 0.;
#endif
            ithread_loc              = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
            Nbre_thread_actif_loc    = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
            E_Int* ipt_ind_dm_thread = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

            E_Float* ipt_CL = iptro_CL[nd];

            if( kimpli == 1 )
              {  
		 E_Int lrhs=1; E_Int lcorner=1; E_Int mjrnewton=1;
                 if (ithread_loc == -1) {continue;}

                 if(lexit_lu == 0 )
                  { 
#include "FastS/Compute/LU/prep_lussor.h"

                    if(barrier_residu==0)
                    {
#include             "FastS/Compute/verrou_lhs.h"
                    }
                    
                    E_Int* ipt_shift_lu   = shift_lu + 6*ntask*Nbre_thread_actif + (ithread-1)*6;

                    E_Float* iptdrodm_out = iptdrodm + shift_zone;

                    if(param_int[nd][LU_MATCH]==1 || param_int[nd][NB_RELAX]>1) iptdrodm_out = ipt_ssortmp_shift;

		    invlu_(nd                     , nitcfg                  , nitrun                ,
			   param_int[nd]          , param_real[nd]          ,
			   ipt_shift_lu           , ipt_ind_dm_thread       , mjrnewton             ,
			   iptrotmp[nd]           , iptro_ssiter[nd]        , iptdrodm + shift_zone , iptdrodm_out,
			   ipti[nd]               , iptj[nd]                , iptk[nd]              ,
			   iptventi[nd]           , iptventj[nd]            , iptventk[nd]          ,
			   iptcoe  + shift_coe    , ipt_ssor_shift          , ssor_size);
                             
                    if(nitrun*nitcfg > 15 and (nitcfg < 3 or nitcfg == nssiter-1)) //protection garbage collector
                      {
#ifdef _OPENMP
                        E_Float lhs_end = omp_get_wtime();
#else
                        E_Float lhs_end = 0.;
#endif
                        E_Int cpu_perzone   =  nssiter*Nbre_thread_actif*2 + nd*(Nbre_thread_actif*2+1);
                        E_Int cells = (ipt_shift_lu[1]-ipt_shift_lu[0]+1)*(ipt_shift_lu[3]-ipt_shift_lu[2]+1)*(ipt_shift_lu[5]-ipt_shift_lu[4]+1);
                        E_Int ith = ithread;
                        timer_omp[ cpu_perzone + 1+ (ith-1)*2 ] +=(lhs_end - lhs_begin)/double(cells);


                      }
                  }//lexit
              } //fin kimpli


            // Selective Frequency Damping
            if(( (nitcfg == nssiter && lssiter_verif==1) || (nitcfg == nssiter-1 && lssiter_verif==0)) && param_int[nd][SFD] == 1)
              {
               sfd_(param_int[nd], param_real[nd], nitrun, ipt_ind_dm_thread, ipt_CL, iptrof[nd], iptcoe + shift_coe, iptvol[nd]);
              }

        }//maillage structure
      }//fin boucle task


