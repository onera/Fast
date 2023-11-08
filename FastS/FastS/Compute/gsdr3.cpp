# include "FastS/fastS.h"
# include "FastC/fastc.h"
# include "FastS/param_solver.h"
# include "connector.h"
# include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

using namespace K_FLD;
using namespace std;

#undef Conservatif
//#define Conservatif
#undef TimeShow
//#define TimeShow

#ifdef TimeShow

#include <iostream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <sstream>

E_Float time_COM=0.0;
E_Float time_init;
#endif

//E_Float time_COM=0.0;
//E_Float time_init;


E_Int K_FASTS::gsdr3( 
  E_Int**& param_int  , E_Float**& param_real ,
  E_Int& nidom        , E_Int& nitrun         , E_Int&  nitcfg    , E_Int&  nitcfg_last, E_Int&  nssiter, E_Int& it_target , E_Int&  first_it,
  E_Int& kimpli       , E_Int& lssiter_verif  , E_Int& lexit_lu   , E_Int& layer_mode  , E_Int& mpi,
  E_Int& nisdom_lu_max, E_Int& mx_nidom       , E_Int& ndimt_flt  ,
  E_Int& threadmax_sdm, E_Int& mx_synchro,
  E_Int& nb_pulse     , 
  E_Float& temps      , E_Float& time_trans,
  E_Int* ipt_ijkv_sdm ,
  E_Int* ipt_ind_dm_omp, E_Int* ipt_topology    , E_Int* ipt_ind_CL    , E_Int* ipt_lok,  E_Int* verrou_lhs,  E_Int& vartype, E_Float* timer_omp,
  E_Int*     iptludic  , E_Int*   iptlumax      ,
  E_Int** ipt_ind_dm   , E_Int** ipt_it_lu_ssdom, E_Int* iptdtloc, 
  E_Float* ipt_VectG   , E_Float* ipt_VectY     , E_Float** iptssor    , E_Float** iptssortmp    , E_Int* ipt_ssor_size, E_Float* ipt_drodmd,
  E_Float* ipt_Hessenberg, E_Float** iptkrylov      , E_Float** iptkrylov_transfer, E_Float* ipt_norm_kry, E_Float** ipt_gmrestmp, E_Float* ipt_givens,
  E_Float*   ipt_cfl,
  E_Float**  iptx            , E_Float**  ipty            , E_Float**    iptz      ,
  E_Float**  iptCellN        , E_Float**  iptCellN_IBC    , E_Int** ipt_degen      ,
  E_Float**& iptro           , E_Float**& iptro_m1        , E_Float**&  iptrotmp   , E_Float**& iptrof       ,
  E_Float**  iptmut          , E_Float*   ipt_mutd        ,
  E_Float**  ipti            , E_Float**  iptj            , E_Float** iptk         , E_Float** iptvol        ,
  E_Float**  ipti0           , E_Float**  iptj0           , E_Float** iptk0        ,
  E_Float**  ipti_df         , E_Float**  iptj_df         , E_Float** iptk_df      , E_Float** iptvol_df     ,
  E_Float**  iptventi        , E_Float**  iptventj        , E_Float** iptventk     ,
  E_Float**& iptrdm          ,
  E_Float*   iptroflt        , E_Float*   iptroflt2       , E_Float*  iptwig       , E_Float*   iptstat_wig  ,
  E_Float*   iptdrodm        , E_Float*   iptcoe          , E_Float*  iptrot       , E_Float**& iptdelta     , E_Float**& iptro_res, E_Float**& iptdrodm_transfer,
  E_Int*&    param_int_tc    , E_Float*& param_real_tc    , E_Int*& linelets_int   , E_Float*& linelets_real , 
  E_Int&     taille_tabs     , E_Float*& stock            , E_Float*& drodmstock   , E_Float*& constk        , E_Float** iptsrc)

 {
#ifdef TimeShow  
  E_Int rank =0; 
  E_Float* ipt_timecount = new E_Float[5];
  ipt_timecount[0:4] = 0.0;
  E_Int nbpointsTot  =0;

  for (E_Int nd = 0; nd < nidom; nd++) { nbpointsTot = nbpointsTot + param_int[nd][ IJKV ]*param_int[nd][ IJKV +1 ]*param_int[nd][ IJKV +2 ]; }
#else
  E_Float* ipt_timecount = NULL;        
#endif

#ifdef TimeShow

#ifdef _MPI
  if(mpi) { MPI_Comm_rank (MPI_COMM_WORLD, &rank); }
#endif

  time_init = 0;
#ifdef _OPENMP  
  time_init = omp_get_wtime();
#endif
    
#endif

      E_Int omp_mode = iptdtloc[8];
      E_Int shift_omp= iptdtloc[11];
      E_Int* ipt_omp = iptdtloc + shift_omp;

      E_Int npass         = 0;
      E_Int ibord_ale     = 1;      // on autorise un calcul optimisee des vitesse entrainement en explicit
      E_Int nptpsi        = 1;
      E_Int balance       = 0;

   //if(nitcfg > 1) return ibord_ale;

      E_Int mx_sszone = mx_nidom/nidom;

      FldArrayI tot( 6*threadmax_sdm); E_Int* ipt_tot  =  tot.begin();



      //
      // choix du tableau de travail en fonction du schema et sous-iteration
      //
      ////  Verifier iptro_CL pour nitcfg> 1 en implicite
      //
      //
      E_Float** iptro_ssiter;
      E_Float** iptro_CL; 

      if (param_int[0][EXPLOC]== 1 and param_int[0][ITYPCP]==2)   //explicit local instationnaire
	   {

	   if (nitcfg%2 != 0){iptro_ssiter = iptro;  iptro_CL = iptrotmp;}
	   else
	     {
	        if (nitcfg != param_int[0][NSSITER]){iptro_ssiter = iptrotmp;  iptro_CL = iptro;}
               else {iptro_ssiter  = iptrotmp; iptro_CL = iptro;}
              }
	   }
      else  // Explicit global ou Implicit
	    {
	       if( nitcfg == 1) { iptro_ssiter = iptro;  iptro_CL = iptrotmp;}
	       else
	       {      
		 if      (param_int[0][ ITYPCP ]  < 2                                                    ) {iptro_ssiter = iptrotmp; iptro_CL = iptrotmp;} // Implicite
		 if      (param_int[0][ ITYPCP ] == 2 && nitcfg%2==0 && nitcfg != param_int[0][ NSSITER ]) {iptro_ssiter = iptrotmp; iptro_CL = iptro;   } // Explicite
		 else if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==0 && nitcfg == param_int[0][ NSSITER ]) {iptro_ssiter = iptrotmp; iptro_CL = iptro;   } // Explicite
		 if      (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg != param_int[0][ NSSITER ]) {iptro_ssiter = iptro;    iptro_CL = iptrotmp;} // Explicite
		 else if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg == param_int[0][ NSSITER ]) {iptro_ssiter = iptro;    iptro_CL = iptrotmp;} // Explicite
	      }
	    }

      E_Int nbtask = ipt_omp[nitcfg-1]; 
      E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

      //mise a jour it_bloc implicit local, pointeur volume en deformable  et  pointeur ssor
#include   "FastS/Compute/ssor.cpp"

//modif Guillaume??
if(nitcfg==1){param_real[0][TEMPS] = 0.0;}

 
/****************************************************
----- Debut zone // omp
****************************************************/

#ifdef _OPENMP
  E_Int Nthread_max  = omp_get_max_threads();
#else
  E_Int Nthread_max  = 1;
#endif

  E_Int iptflux      = param_int[0][IBC_PT_FLUX];
  E_Int Nfamily      = 0;
  E_Int neqFlu       = 26;
  if(iptflux != -1) Nfamily  = param_int[0][iptflux];

  E_Float masse[3];
  E_Float debit[Nthread_max*2*6];
  E_Float flux[Nthread_max*neqFlu*Nfamily];
  E_Float ro_corr; E_Float varMasse;

#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
    E_Float rhs_begin        = omp_get_wtime();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
    E_Float rhs_begin       = 0;
#endif
   //E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;

   E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
   E_Int socket          = (ithread-1)/thread_parsock +1;
   E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

   E_Int* ipt_topology_socket    = ipt_topology       + (ithread-1)*3;
   E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm       + (ithread-1)*3;
   E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
   E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
   E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
   E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;
   E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;
   E_Int* ipt_ind_dm_omp_thread  = ipt_ind_dm_socket  + 6;

   E_Int* ipt_nidom_loc, nb_subzone;
     /****************************************************
      -----Boucle sous-iteration
     ****************************************************/
        E_Int nbtask = ipt_omp[nitcfg-1]; 
        E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

        if( nitcfg == 1)
                  {
                     //mise a jour metric et vent ale zone cart et 3dhom(3dfull et 2d a la volee)
                     for (E_Int ntask = 0; ntask < nbtask; ntask++)
                         {
                            E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
                            E_Int nd = ipt_omp[ pttask ];
                           if(param_int[nd][LALE]==1) //maillage indeformable
                             {
                               mjr_ale_3dhomocart_(nd, param_int[nd] ,   param_real[nd]   ,
                                                  socket             ,  Nbre_socket       , ithread_sock        , thread_parsock,
                                                  ipt_ind_dm_socket  , ipt_topology_socket,
                                                  iptx[nd]           , ipty[nd]           , iptz[nd]            ,
                                                  ipti[nd]           , iptj[nd]           , iptk[nd]            ,
                                                  ipti0[nd]          , iptj0[nd]          , iptk0[nd]           , iptvol[nd] ,
                                                  iptventi[nd]       , iptventj[nd]       , iptventk[nd]        );
                                 //modifier mjr_ale_3dhomocart_ pour faire sauter barrier
                                 #pragma omp barrier
                             }
                         }//zone


                   // calcul metric si maillage deformable
                   //  
#include           "FastS/Metric/cp_metric.cpp"
                   }


        //---------------------------------------------------------------------
        // -----Boucle sur num.les domaines de la configuration
        // ---------------------------------------------------------------------
        E_Int shift_zone=0; E_Int shift_wig=0; E_Int shift_coe=0; E_Int shift_mu=0; E_Int nd_current=0;
        E_Float rhs_end=0; E_Int ncells=0;

        if (param_int[0][IMPLICITSOLVER] == 1 && layer_mode == 1) { ipt_norm_kry[ithread-1]=0.; }


        for (E_Int ntask = 0; ntask < nbtask; ntask++)
          {
             E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
             E_Int nd = ipt_omp[ pttask ];

             E_Int lmin = 10;
             if (param_int[nd][ITYPCP] == 2) lmin = 4;

             shift_zone=0; shift_wig=0; shift_coe=0;
             for (E_Int n = 0; n < nd; n++)
             {
              shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
              shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
              if(param_int[n][ KFLUDOM ]==2){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*3;}
              if(param_int[n][ SLOPE   ]==6){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*3;}
              if(param_int[n][ SLOPE   ]==7){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*3;}
             }
#include "FastS/Compute/rhs.cpp"
          } 
#ifdef Conservatif
#include   "FastS/Compute/cp_debitIBM.cpp"
#endif

#ifdef _WIN32
#pragma omp barrier
#endif
          //
          //timer pour omp "dynamique"
          //
#ifdef _OPENMP  
          rhs_end = omp_get_wtime();
#endif
          E_Int cpu_perthread = (ithread-1)*2 +(nitcfg-1)*Nbre_thread_actif*2;
          timer_omp[ cpu_perthread] += rhs_end - rhs_begin;
          //if(ncells !=0 )printf(" time rhs= %g %d  %d %d %d  \n",(rhs_end - rhs_begin)/float(ncells) , ithread, ncells, nitcfg, nssiter);
	  //if(ithread==1) printf("========= \n");
          //if(ithread==1) printf(" time rhs= %g %g  %d %d  \n",timer_omp[ cpu_perthread], rhs_end - rhs_begin , cpu_perthread, nitcfg);

          if (param_int[0][IMPLICITSOLVER] == 1 && layer_mode == 1)
	  {
#include   "FastS/Compute/Linear_solver/lhs.cpp"
          }


          // LUSGS
          else
	  { 
#include   "FastS/Compute/lhs.cpp"

#ifdef Conservatif
//#include   "FastS/Compute/conservatif.cpp"
#endif
          }


          //
          //finalisation timer pour omp "dynamique"
          //
#ifdef _OPENMP  
          E_Float     lhs_end = omp_get_wtime();
#else  
          E_Float     lhs_end = 0.;
#endif
          timer_omp[ cpu_perthread +1 ] += lhs_end- rhs_end;

} // Fin zone // omp


#ifdef TimeShow
#ifdef _OPENMP  
     time_COM = omp_get_wtime();
     time_init = omp_get_wtime();
#endif 
#endif 


#ifdef _OPENMP  
     E_Float     trans_deb = omp_get_wtime();
#else  
     E_Float     trans_deb = 0.;
#endif

  //
  //
  //FillGhostcell si mise a jour necessaire et transfer dans C layer 
  //

E_Int autorisation_bc[nidom];
for (E_Int nd = 0; nd < nidom; nd++) { autorisation_bc[nd]=1;}

E_Int nitcfg_stk = nitcfg;
if(lexit_lu ==0 && layer_mode==1)
{   
  //remplissage ghost transfert
  for (E_Int ipass = 0; ipass < 1; ipass++)
  {
  #include "FastS/Compute/transfert_multiblock.cpp"
  }
  E_Int cycl;
 
  for (E_Int nd = 0; nd < nidom; nd++)
    {
     // flag pour transfert optimiser explicit local
     autorisation_bc[nd]=0;
     if (param_int[0][EXPLOC] == 1)    //if (rk == 3 and exploc == 2) 
      {
	cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];
	 
        // modif de nitcfg pour appliquer les BC partout 
	if      (nitcfg_stk%cycl == cycl/2 -1       and cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}
	else if (nitcfg_stk%cycl == cycl/2 + cycl/4 and cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}	    
	else if (nitcfg_stk%cycl == cycl-1          and cycl != 4 ){ nitcfg = 1; autorisation_bc[nd] = 1;}
	else if((nitcfg_stk%cycl == 1 or nitcfg_stk%cycl == cycl/2  or nitcfg_stk%cycl== cycl-1) and cycl == 4 ) { nitcfg = 1; autorisation_bc[nd] = 1; }
      } 
     else {autorisation_bc[nd] = 1;}
  }//loop zone
} //lexit  

#ifdef _OPENMP  
    time_trans = time_trans + omp_get_wtime()-trans_deb;
#else  
    time_trans = 0.;
#endif

//calcul sequentiel flux IBM et data CLs (lund, wallmodel) a toutes les ssiter
for (E_Int nd = 0; nd < nidom; nd++)
  {
     #include  "FastS/Compute/data_ibm_bc.cpp"
  }//loop zone
 
//filtre domaine a traiter en IBM pour optimiser loop // for
/*E_Int nd_ibc[nidom]; E_Int nidom_ibc =0;
for (E_Int nd = 0; nd < nidom; nd++)
  {
   E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
   E_Int  nb_subzone    = ipt_nidom_loc [nitcfg_stk-1];                                       //nbre sous-zone a la sousiter courante

   //if (autorisation_bc[nd] == 1 && nb_subzone > 0)
   // {
      E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
      if (pt_flu != -1) { nd_ibc[nidom_ibc]=nd; nidom_ibc +=1; }
   // }
  }
*/
#ifdef _OPENMP  
     trans_deb = omp_get_wtime();
#else 
     trans_deb = 0.;
#endif

E_Int lrhs=0; E_Int lcorner=0; 
#pragma omp parallel default(shared)
 {
#ifdef _OPENMP
   E_Int  ithread           = omp_get_thread_num() +1;
   E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
   E_Int ithread = 1;
   E_Int Nbre_thread_actif = 1;
#endif

   E_Int Nbre_thread_actif_loc, ithread_loc;
   //
   //Apply BC (parcour Zones) + reinitialisation verrou pour calcul rhs
   //
   if(lexit_lu ==0 && layer_mode==1)
   { 
    #include   "FastS/Compute/bcs.cpp"
   }

 }//fin zone omp

#ifdef _OPENMP  
     time_trans = time_trans + omp_get_wtime()-trans_deb;
#else 
     time_trans = 0.;
#endif

  nitcfg = nitcfg_stk;

//#ifdef Conservatif
//#include "FastS/Compute/postIBM.cpp"
//#endif

    //
    //omp "dynamic" balance
    //
//    if(nitrun%5==0)
//    {  
//#include "FastC/HPC_LAYER/OPTIMIZE_DISTRIB_OMP.h"
//
//    }

   delete [] ipt_timecount;
   return ibord_ale;

}
