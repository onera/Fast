# include "FastS/fastS.h"
# include "fastc.h"
# include "FastS/param_solver.h"
# include "connector.h"
# include <string.h>
//# include <CMP/include/pending_message_container.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

using namespace K_FLD;
using namespace std;
using namespace K_CONNECTOR;
using namespace K_FASTC;

#undef TimeShow
// #define TimeShow

#ifdef TimeShow

#include <iostream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <sstream>

E_Float time_COM=0.0;
E_Float time_init;
#endif

 void K_FASTS::matvec( 
  E_Int**& param_int  , E_Float**& param_real , E_Int& no_vect_test,
  E_Int& nidom        , E_Int& nitrun         , E_Int&  nitcfg    , E_Int&  nssiter, E_Int& it_target , E_Int&  first_it,
  E_Int& kimpli       , E_Int& lssiter_verif  , E_Int& lexit_lu   , E_Int& omp_mode, E_Int& layer_mode, E_Int& mpi,
  E_Int& nisdom_lu_max, E_Int& mx_nidom       , E_Int& ndimt_flt  ,
  E_Int& threadmax_sdm, E_Int& mx_synchro,
  E_Int& nb_pulse     , 
  E_Float& temps,
  E_Int* ipt_ijkv_sdm  ,
  E_Int* ipt_ind_dm_omp, E_Int* ipt_topology    , E_Int* ipt_ind_CL    , E_Int* ipt_lok,  E_Int* verrou_lhs,  E_Int& vartype, E_Float* timer_omp,
  E_Int*     iptludic  , E_Int*   iptlumax      , E_Int* ipt_omp       ,
  E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom  ,
  E_Float* ipt_VectG      , E_Float* ipt_VectY   , E_Float** iptssor       , E_Float** iptssortmp    , E_Int* ipt_ssor_size, E_Float* ipt_drodmd,
  E_Float* ipt_Hessenberg, E_Float** iptkrylov      , E_Float** iptkrylov_transfer, E_Float* ipt_norm_kry, E_Float** ipt_gmrestmp, E_Float* ipt_givens,
  E_Float*   ipt_cfl,
  E_Float**  iptx            , E_Float**  ipty            , E_Float**    iptz      ,
  E_Float**  iptCellN        , E_Float**  iptCellN_IBC    ,
  E_Float**& iptro           , E_Float**& iptro_m1        , E_Float**&  iptrotmp   , E_Float**& iptrof       ,
  E_Float**  iptmut          , E_Float*   ipt_mutd        ,
  E_Float**  ipti            , E_Float**  iptj            , E_Float** iptk         , E_Float** iptvol        ,
  E_Float**  ipti0           , E_Float**  iptj0           , E_Float** iptk0        ,
  E_Float**  ipti_df         , E_Float**  iptj_df         , E_Float** iptk_df      , E_Float** iptvol_df     ,
  E_Float**  iptventi        , E_Float**  iptventj        , E_Float** iptventk     ,
  E_Float**& iptrdm          ,
  E_Float*   iptroflt        , E_Float*   iptroflt2       , E_Float*  iptwig       , E_Float*   iptstat_wig  ,
  E_Float*   iptdrodm        , E_Float*   iptcoe          , E_Float*  iptrot       , E_Float**& iptdelta , E_Float**& iptro_res,
  E_Float**& iptdrodm_transfer, 
  E_Int*&   param_int_tc , E_Float*& param_real_tc, E_Int*& linelets_int, E_Float*& linelets_real)

 {




      E_Int npass         = 0;
      E_Int ibord_ale     = 1;      // on autorise un calcul optimisee des vitesse entrainement en explicit
      E_Int nptpsi        = 1;
      E_Int balance       = 0;

      FldArrayI tot( 6*threadmax_sdm); E_Int* ipt_tot  =  tot.begin();


      //
      // choix du tableau de travail en fonction du schema et sous-iteration
      //
      ////  Verifier iptro_CL pour nitcfg> 1 en implicite
      //
      //
      E_Float** iptro_ssiter;
      E_Int ishift,lfwmean;
      E_Float** iptro_CL; 
      if (param_int[0][EXPLOC]== 1 and param_int[0][ITYPCP]==2) // Explicit local
      {
         if (nitcfg%2 != 0){iptro_ssiter = iptro;  iptro_CL = iptrotmp; ishift  =1; lfwmean  = 1; }
	 else
	     {
	        if (nitcfg != param_int[0][NSSITER]){iptro_ssiter = iptrotmp;  iptro_CL = iptro; ishift  = -1; lfwmean = 0;}
                else {iptro_ssiter  = iptrotmp; iptro_CL = iptro;}
             }
      }
      else  // Explicit global ou Implicit
      {

	 if( nitcfg == 1) { iptro_ssiter = iptro;  iptro_CL = iptrotmp; ishift  = 1; lfwmean  = 1; }
	 else
	      {           iptro_ssiter = iptrotmp;  iptro_CL = iptro; ishift =-1; lfwmean  = 0;
                          if (param_int[0][ ITYPCP ] == 2 && nitcfg == 3) {iptro_ssiter  = iptro; iptro_CL = iptrotmp;}
			  if (param_int[0][ ITYPCP ] < 2                ) {                       iptro_CL = iptrotmp;}
	      }
      }

      for (E_Int nd = 0; nd < nidom; nd++)
	  {
	   E_Float* matvec_x = iptkrylov[nd] +   no_vect_test   *param_int[nd][NEQ] * param_int[nd][NDIMDX];
           E_Float* matvec_b = iptkrylov[nd] +  (no_vect_test+1)*param_int[nd][NEQ] * param_int[nd][NDIMDX];

	   iptkrylov_transfer[nd] = matvec_x;
          }

      E_Int rk = 3;
      E_Int exploc = 0;
      E_Int numpassage = 1;
      //Raccord X
      E_Float* ipt_timecount = NULL;
if( param_int_tc != NULL)
  { 
      K_FASTC::setInterpTransfersFast(iptkrylov_transfer, vartype, param_int_tc,
			      param_real_tc, param_int, param_real, ipt_omp, linelets_int, linelets_real,
         		      it_target, nidom, ipt_timecount, mpi, nitcfg , nssiter, rk, exploc, numpassage);
  } 
   


/****************************************************
----- Debut zone // omp
****************************************************/
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
   //E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;
   if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
   else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

   E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
   E_Int socket          = (ithread-1)/thread_parsock +1;
   E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

   E_Int nbtask; E_Int ptiter;
   E_Int shift_zone=0; E_Int shift_wig=0; E_Int shift_coe=0; E_Int shift_mu=0;
   E_Int lrhs = 2; E_Int lcorner = 0; E_Int npass = 0;
    E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
    E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
    E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
    E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;
#include  "HPC_LAYER/OMP_MODE_BEGIN.h"

           E_Int ierr = BCzone_d(nd, lrhs , nitcfg, lcorner, param_int[nd], param_real[nd], npass,
	           		 ipt_ind_dm_loc, ipt_ind_dm_thread, 
			         ipt_ind_CL_thread, ipt_ind_CL119 ,  ipt_ind_CLgmres, ipt_shift_lu,
			         iptro_ssiter[nd],
		    	         ipti[nd]     , iptj[nd]            , iptk[nd],
			         iptx[nd]     , ipty[nd]            , iptz[nd],
			         iptventi[nd] , iptventj[nd]        , iptventk[nd],
			         iptkrylov_transfer[nd] );

            correct_coins_(nd,  param_int[nd], ipt_ind_dm_thread , iptkrylov_transfer[nd]);

#include    "HPC_LAYER/OMP_MODE_END.h"


#pragma omp barrier

    //
    //
    // produit matrice:vecteur.
    // In: rop_ssiter, rop_ssiter_d
    // Out: drodmd
    ipt_norm_kry[ithread-1]=0;

    E_Int* ipt_topology_socket    = ipt_topology       + (ithread-1)*3;
    E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm       + (ithread-1)*3;
    E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;
    E_Int* ipt_ind_dm_omp_thread  = ipt_ind_dm_socket  + 6;

#include  "HPC_LAYER/OMP_MODE_BEGIN.h"

       E_Float* krylov_in    = NULL;
       E_Float* rop_ssiter_d = iptkrylov[nd] +  no_vect_test    * param_int[nd][NEQ] * param_int[nd][NDIMDX]; //matvec_x
       E_Float* ipt_drodmd   = iptkrylov[nd] + (no_vect_test+1) * param_int[nd][NEQ] * param_int[nd][NDIMDX]; //matvec_out

       E_Int mjr_dt =1;
#include "Compute/Linear_solver/dRdp_tapenade.cpp"

	//shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
	// Warning
	// Warning a blinder
	// Warning
	//shift_zone = 0;
#include    "HPC_LAYER/OMP_MODE_END.h"

  }// omp


  for (E_Int nd = 0; nd < nidom; nd++)
  {
   E_Float* matvec_b = iptkrylov[nd] +  (no_vect_test+1)*param_int[nd][NEQ] * param_int[nd][NDIMDX];

   iptkrylov_transfer[nd] = matvec_b;
  }


      //Raccord X
if( param_int_tc != NULL)
  { 
      K_FASTC::setInterpTransfersFast(iptkrylov_transfer, vartype, param_int_tc,
			      param_real_tc, param_int, param_real, ipt_omp, linelets_int, linelets_real,
         		      it_target, nidom, ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage);
  } 
   

 }
