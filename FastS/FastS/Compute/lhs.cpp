if( kimpli == 1  && param_int[0][LU_MATCH]==1)
  { 
   #pragma omp master
    { //Raccord V0

      setInterpTransfersFastS(iptdrodm_transfer, ndimdx_transfer, param_int_tc,
                           param_real_tc, param_int_ibc, param_real_ibc, linelets_int, linelets_real,param_real[0][PRANDT],
                           it_target, nidom, ipt_timecount, mpi);
    }
    #pragma omp barrier
  }


            shift_zone=0; shift_coe=0; nd_current=0;
            E_Float* ipt_ssor_shift; E_Float* ipt_ssortmp_shift; E_Int ssor_size;
            for (E_Int nd = 0; nd < nidom; nd++)
            {
               E_Float* ipt_CL = iptro_CL[nd];

#include       "HPC_LAYER/OMP_MODE_BEGIN.h"
                      //
                      //verrou rhs
                      //
                      E_Int type = 2;
                      for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
                      { 
                           E_Int* verrou_lhs_thread= verrou_lhs + nd_current*Nbre_thread_actif +th;
                           verrou_c_( verrou_lhs_thread, type );
                      }

                      //#pragma omp flush
                      E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX];
                      //flush_real_( size , iptdrodm + shift_zone);
                      if(nitcfg==1)
                      {
                        size = param_int[nd][NEQ_COE]*param_int[nd][NDIMDX];
                        flush_real_( size , iptcoe + shift_coe);
                      }
                      //size       = param_int[nd][NDIMDX];
                      //flush_real_( size , iptmut[nd]);

                     //sortie de la carte d residu du Newton
#include             "Compute/residus_navier.h"

                     if( kimpli == 1 )
                      {  //
                         // CL sur rhs pour implicitation
			 E_Int lrhs=1; E_Int lcorner=1; E_Int mjrnewton=1;
                         BCzone( nd, lrhs, lcorner,
                                 param_int[nd], param_real[nd],
                                 npass,
                                 ipt_ind_dm_loc         , ipt_ind_dm_thread      ,
                                 ipt_ind_CL_thread      , ipt_ind_CL119          , ipt_ind_CLgmres, ipt_shift_lu,
                                 iptdrodm + shift_zone  , ipti[nd]               , iptj[nd]       , iptk[nd]          ,
                                 iptx[nd]               , ipty[nd]               , iptz[nd]       ,
                                 iptventi[nd]           , iptventj[nd]           , iptventk[nd]   , iptrotmp[nd]);
                 
                         if(lcorner  == 0 )correct_coins_(nd, param_int[nd], ipt_shift_lu , iptdrodm + shift_zone );

                         if(lssiter_verif ==1  && nd_subzone ==0 && omp_mode==1 && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
                         {
                          E_Int type = 2;
                          for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
                            { 
                              E_Int* verrou_lhs_thread= verrou_lhs + (mx_nidom + nd_current)*Nbre_thread_actif +th;
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
			   }
                      } //fin kimpli


                     // Selective Frequency Damping
                     if(( (nitcfg == nssiter && lssiter_verif==1) || (nitcfg == nssiter-1 && lssiter_verif==0)) && param_int[nd][SFD] == 1)
                     {
                       sfd_(param_int[nd], param_real[nd], nitrun, ipt_ind_dm_thread, ipt_CL, iptrof[nd], iptcoe + shift_coe, iptvol[nd]);
                     }

                   nd_current +=1;
                     
#include       "HPC_LAYER/OMP_MODE_END.h"
                 //  } //fin boucle sous-zone

             shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
             shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
            }//fin boucle zone


