if( kimpli == 1  && param_int[0][LU_MATCH]==1 && param_int_tc != NULL)
  { 
   #pragma omp master
    { //Raccord V0
      //setInterpTransfersFastS(iptdrodm_transfer, vartype, param_int_tc,
      //                   param_real_tc, param_int, param_real, linelets_int, linelets_real,
      // 			 it_target, nidom, ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage);
      
      E_Int rk           = param_int[0][RK];
      E_Int exploc       = param_int[0][EXPLOC];
      E_Int numpassage   = 1;

      K_FASTC::setInterpTransfersFast(iptdrodm_transfer, vartype, param_int_tc,
                         param_real_tc, param_int, param_real, linelets_int, linelets_real,
      			 it_target, nidom, ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage);
    }
    #pragma omp barrier
  }


            shift_zone=0; shift_coe=0; nd_current=0;
            E_Float* ipt_ssor_shift; E_Float* ipt_ssortmp_shift; E_Int ssor_size;
            for (E_Int nd = 0; nd < nidom; nd++)
            {

             if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
             {
               E_Float* ipt_CL = iptro_CL[nd];

#ifdef _OPENMP
               E_Float lhs_begin = omp_get_wtime();
#else
               E_Float lhs_begin = 0.;
#endif
#include       "HPC_LAYER/OMP_MODE_BEGIN.h"
                      //
                      //verrou rhs
                      //
                      E_Int type = 2;
                      //si calcul residu on verifie que toutes les souszones sont pretes
                      if(lssiter_verif ==1  && nd_subzone ==0 && omp_mode==1 && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
                      {
                        if (nb_subzone !=1)// pour 1 drodm OK pour calcul residu
                        {
                          E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX];
                          flush_real_( size , iptdrodm + shift_zone);
                        }
                        for (E_Int nd_subzone_loc = 0; nd_subzone_loc < nb_subzone; nd_subzone_loc++)
                        {
                           E_Int       Ptomp = param_int[nd][PT_OMP];
                           E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
                           E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone_loc];

                           E_Int Nbre_thread_actif_loc1 = param_int[nd][ PtZoneomp  + Nbre_thread_actif ];
                           E_Int ithread_loc1           = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
  
                           if (ithread_loc1 == -1) { continue;}

                           for (E_Int th = 0; th < Nbre_thread_actif_loc1; th++) { 
                             E_Int* verrou_lhs_thread= verrou_lhs + (nd_current+nd_subzone_loc)*Nbre_thread_actif +th;
                             verrou_c_( verrou_lhs_thread, type );
                           }
                        }
                      }
                      else
                      { //on verifie uniquememnt la souszone courante: nd_current
                         for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) { 
                           E_Int* verrou_lhs_thread= verrou_lhs + nd_current*Nbre_thread_actif +th;
                           verrou_c_( verrou_lhs_thread, type );
                         }
                      }
                      if(nitcfg==1)
                      {
                        E_Int size = param_int[nd][NEQ_COE]*param_int[nd][NDIMDX];
                        flush_real_( size , iptcoe + shift_coe);
                      }

                     //sortie de la carte d residu du Newton
#include             "Compute/residus_navier.h"

                     if( kimpli == 1 )
                      {  //
                         // CL sur rhs pour implicitation
			 E_Int lrhs=1; E_Int lcorner=1; E_Int mjrnewton=1;
                         K_FASTS::BCzone( nd, lrhs, nitcfg, lcorner, 
                                          param_int[nd], param_real[nd],
                                          npass,
                                          ipt_ind_dm_loc         , ipt_ind_dm_thread      ,
                                          ipt_ind_CL_thread      , ipt_ind_CL119          , ipt_ind_CLgmres, ipt_shift_lu,
                                          iptdrodm + shift_zone  , ipti[nd]               , iptj[nd]       , iptk[nd]          ,
                                          iptx[nd]               , ipty[nd]               , iptz[nd]       ,
                                          iptventi[nd]           , iptventj[nd]           , iptventk[nd]   , iptrotmp[nd], iptmut[nd]);
                 
                         if(lcorner  == 0 )correct_coins_(nd, param_int[nd], ipt_shift_lu , iptdrodm + shift_zone );

                         if(lssiter_verif ==1  && nd_subzone ==0 && omp_mode==1 && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
                         {
                          E_Int type = 2;
                          for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
                            { 
                              E_Int* verrou_lhs_thread= verrou_lhs + (mx_nidom + nd_current)*Nbre_thread_actif +th;
                          //if(nd==1 and nitrun==881 && nitcfg==1 ) printf("HUM_WAIT %d %d %d %d th= %d\n", (mx_nidom+nd_current)*Nbre_thread_actif, ithread, ithread_loc, Nbre_thread_actif_loc, th);
                              verrou_c_( verrou_lhs_thread, type );
                            }
                         } //sinon residu pas bon en omp_mode=1

                         if(lexit_lu == 0 )
			   { 
#include "Compute/LU/prep_lussor.h"

                             E_Float* iptdrodm_out = iptdrodm + shift_zone;
                             if(param_int[nd][LU_MATCH]==1 || param_int[nd][NB_RELAX]>1) iptdrodm_out = ipt_ssortmp_shift;

			     invlu_(nd                     , nitcfg                  , nitrun                ,
				    param_int[nd]          , param_real[nd]          ,
				    ipt_shift_lu           , ipt_ind_dm_thread       , mjrnewton             ,
				    iptrotmp[nd]           , iptro_ssiter[nd]        , iptdrodm + shift_zone , iptdrodm_out,
				    ipti[nd]               , iptj[nd]                , iptk[nd]              ,
				    iptventi[nd]           , iptventj[nd]            , iptventk[nd]          ,
				    iptcoe  + shift_coe    , ipt_ssor_shift          , ssor_size);
                             
                             if(nitrun*nitcfg > 15) //protection garbage collector
                             {
#ifdef _OPENMP
                               E_Float lhs_end = omp_get_wtime();
#else
                               E_Float lhs_end = 0.;
#endif
                               E_Int cpu_perzone   =  nssiter*Nbre_thread_actif*2 + nd*(Nbre_thread_actif+1);
                               E_Int cells = (ipt_shift_lu[1]-ipt_shift_lu[0]+1)*(ipt_shift_lu[3]-ipt_shift_lu[2]+1)*(ipt_shift_lu[5]-ipt_shift_lu[4]+1);
                               E_Int ith = ithread;
                               if (omp_mode == 1) ith = ithread_loc;
                               timer_omp[ cpu_perzone + ith ] +=(lhs_end - lhs_begin)/double(cells);
                             }
                             //if(ithread==1)printf("cpu1= %g \n",(lhs_end - lhs_begin)/double(cells) );
			   }
                      } //fin kimpli


                     // Selective Frequency Damping
                     if(( (nitcfg == nssiter && lssiter_verif==1) || (nitcfg == nssiter-1 && lssiter_verif==0)) && param_int[nd][SFD] == 1)
                     {
                       sfd_(param_int[nd], param_real[nd], nitrun, ipt_ind_dm_thread, ipt_CL, iptrof[nd], iptcoe + shift_coe, iptvol[nd]);
                     }

                   nd_current +=1;
                     
#include       "HPC_LAYER/OMP_MODE_END.h"

             }//maillage structure

             shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
             shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
            }//fin boucle zone


