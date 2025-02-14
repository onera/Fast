# include "FastS/fastS.h"
//# include "fastP.h"
//# include "FastLBM/fastLBM.h"
//# include "FastASLBM/fastLBM.h"
# include "FastC/fastc.h"
# include "Fast/fast.h"
# include "Fast/param_solver.h"
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
#undef TimeShow

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

E_Int K_FAST::gsdr3(E_Int**& param_int          , E_Float**& param_real       , E_Int& nidom                , E_Int& nitrun           ,
		    E_Int&  nitcfg              , E_Int&  nitcfg_last         , E_Int&  nssiter             , E_Int& it_target        ,
		    E_Int&  first_it            , E_Int& kimpli               , E_Int& lssiter_verif        , E_Int& lexit_lu         ,
		    E_Int& layer_mode           , E_Int& mpi                  , E_Int& nisdom_lu_max        ,
		    E_Int& mx_nidom             , E_Int& ndimt_flt            , E_Int& threadmax_sdm        , E_Int& mx_synchro       ,
		    E_Int& nb_pulse             , E_Float& temps              , E_Int* ipt_ijkv_sdm         , E_Int* ipt_ind_dm_omp   ,  E_Int* iptdtloc    ,
		    E_Int* ipt_topology         , E_Int* ipt_ind_CL           , E_Int* shift_lu             , E_Int* ipt_lok          , E_Int* verrou_lhs       ,
		    E_Int& vartype              , E_Float* timer_omp          , E_Int*     iptludic         , E_Int*   iptlumax       ,
		    E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom     , E_Int** ipt_ng_pe           , E_Int** ipt_nfconn      ,
		    E_Int** ipt_nfindex         , E_Float* ipt_VectG          , E_Float* ipt_VectY          , E_Float** iptssor       ,
		    E_Float** iptssortmp        , E_Int* ipt_ssor_size        , E_Float* ipt_drodmd         , E_Float* ipt_Hessenberg ,
		    E_Float** iptkrylov         , E_Float** iptkrylov_transfer, E_Float* ipt_norm_kry       , E_Float** ipt_gmrestmp  ,
		    E_Float* ipt_givens         , E_Float*   ipt_cfl          , E_Float**  iptx             , E_Float**  ipty         ,
		    E_Float**    iptz           , E_Float**  iptCellN         , E_Float**  iptCellN_IBC     , E_Float** iptFltrN      , E_Float** iptSpongeCoef, E_Int** ipt_degen   ,
		    E_Float**& iptro            , E_Float**& iptro_m1         , E_Float**&  iptrotmp        , E_Float**& iptrof       , E_Float**& iptS     ,  E_Float**& iptPsiG,
		    E_Float**  iptmut           , E_Float*   ipt_mutd         , E_Float**  ipti             , E_Float**  iptj         ,
		    E_Float** iptk              , E_Float** iptvol            , E_Float**  ipti0            , E_Float**  iptj0        ,
		    E_Float** iptk0             , E_Float**  ipti_df          , E_Float**  iptj_df          , E_Float** iptk_df       ,
		    E_Float** iptvol_df         , E_Float**  iptventi         , E_Float**  iptventj         , E_Float** iptventk      ,
		    E_Float**& iptrdm           , E_Float*   iptroflt         , E_Float*   iptroflt2        , E_Float*  iptgrad       ,
		    E_Float*  iptwig            , E_Float*   iptstat_wig      , E_Float* iptflu             , E_Float*   iptdrodm     ,
		    E_Float*   iptcoe           , E_Float*  iptrot            , E_Float**& iptdelta         , E_Float**& iptro_res    ,
		    E_Float**& iptdrodm_transfer, E_Int*&    param_int_tc     , E_Float*& param_real_tc     , E_Int*& linelets_int    ,
		    E_Float*& linelets_real     , E_Int&     taille_tabs      , E_Float*& stock             , E_Float*& drodmstock    ,
		    E_Float*& constk            , E_Float** iptsrc            , E_Float* f_horseq           , E_Float* a1_pr          ,
            E_Float* a1_fd              , E_Float* a1_hrr             , E_Float* aneq_o3            , E_Float* psi_corr       ,  E_Int& flag_NSLBM)
{
  E_Float*  feq          = iptdrodm;
  E_Float** feq_transfer = iptdrodm_transfer;

#ifdef TimeShow
  E_Float* ipt_timecount = new E_Float[5];
  ipt_timecount[0:4] = 0.0;
  E_Int nbpointsTot =0;

  for (E_Int nd = 0; nd < nidom; nd++){ nbpointsTot = nbpointsTot + param_int[nd][ IJKV ]*param_int[nd][ IJKV +1 ]*param_int[nd][ IJKV +2 ]; }
#else
  E_Float* ipt_timecount = NULL;
#endif

E_Int rank =0;
#ifdef TimeShow
#ifdef _MPI
  if(mpi){ MPI_Comm_rank (MPI_COMM_WORLD, &rank);  }
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

  E_Int mx_sszone = mx_nidom/nidom;

  FldArrayI tot( 6*threadmax_sdm); E_Int* ipt_tot  =  tot.begin();


  //
  // choix du tableau de travail en fonction du schema et sous-iteration
  //
  ////  Verifier iptro_CL pour nitcfg> 1 en implicite
  //
  //

  E_Int ishift,lfwmean;
  E_Float** iptro_ssiter= new E_Float*[nidom];
  E_Float** iptro_CL    = new E_Float*[nidom*4];
  //E_Float** iptro_CL    = new E_Float*[nidom];
  E_Float** iptrom_CL   = new E_Float*[nidom];
  E_Float** iptroS_CL   = new E_Float*[nidom];
  E_Float** iptpsiG_CL;

  E_Int rk     =  param_int[0][RK];
  E_Int exploc = param_int[0][EXPLOC];
  E_Int numpassage = 1;

  for (E_Int nd = 0; nd < nidom; nd++){
    //NS ou LBM
    if ( param_int[nd][IFLOW]!= 4){
      // explicit local
      if(param_int[0][EXPLOC]== 1 && param_int[0][ITYPCP]==2){
	if (nitcfg%2 != 0){iptro_ssiter[nd] = iptro[nd];  iptro_CL[nd] = iptrotmp[nd]; ishift  =1; lfwmean  = 1; }
	else {
	  if (nitcfg != param_int[0][NSSITER]){iptro_ssiter[nd] = iptrotmp[nd]; iptro_CL[nd] = iptro[nd]; ishift  = -1; lfwmean = 0;}
	  else                                {iptro_ssiter[nd] = iptrotmp[nd]; iptro_CL[nd] = iptro[nd];}
	}

      }
      else  { // Explicit global ou Implicit
	if( nitcfg == 1) { iptro_ssiter[nd] = iptro[nd];  iptro_CL[nd] = iptrotmp[nd]; ishift  = 1; lfwmean  = 1; }
	else{
	  if (param_int[0][ ITYPCP ] < 2                )                                             {iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptrotmp[nd]; ishift =-1; lfwmean  = 0;} // Implicit
	  if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==0 && nitcfg != param_int[nd][ NSSITER ] )      {iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptro[nd];    ishift =-1; lfwmean  = 0;} // Explicite
	  else if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==0 && nitcfg == param_int[nd][ NSSITER ] ) {iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptro[nd];                             } // Explicite
	  if      (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg != param_int[nd][ NSSITER ])  {iptro_ssiter[nd] = iptro[nd];     iptro_CL[nd] = iptrotmp[nd]; ishift =-1; lfwmean  = 0;} // Explicite
	  else if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg == param_int[nd][ NSSITER ])  {iptro_ssiter[nd] = iptro[nd];     iptro_CL[nd] = iptrotmp[nd];                          } // Explicite
	}
      }
    }
    else{ // param_int[nd][IFLOW]== 4
      iptro_ssiter[nd]      = iptrotmp[nd];
      iptro_CL[nd]          = iptro[nd];
      iptro_CL[nd+   nidom] = iptrotmp[nd];
      iptro_CL[nd+ 2*nidom] = iptS[nd];
      iptro_CL[nd+ 3*nidom] = iptPsiG[nd];
      iptrom_CL[nd]      = iptro[nd];
      iptroS_CL[nd]      = iptS[nd];
      iptpsiG_CL         = iptPsiG;
    }
  }//zone

  //
  //Calcul taille tableau ssor par thread et mise a jour Nombre sous_iter pour implicit
  //
  E_Int nbtask = ipt_omp[nitcfg-1];
  E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];
#include   "../FastS/FastS/Compute/ssor.cpp"

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
  if (iptflux != -1) Nfamily  = param_int[0][iptflux];

#ifdef _Conservatif
  E_Float masse[3];
  E_Float debit[Nthread_max*2*6];
  E_Float flux[Nthread_max*neqFlu*Nfamily];
  E_Float ro_corr; E_Float varMasse;
#endif

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
   if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
   else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

   E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
   E_Int socket          = (ithread-1)/thread_parsock +1;
   E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

   E_Int* ipt_topology_socket    = ipt_topology       + (ithread-1)*3;
   E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm       + (ithread-1)*3;
   E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
   E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
   E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
   E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;
   E_Int* ipt_ind_dm_omp_thread  = ipt_ind_dm_socket  + 6;

   E_Int* ipt_nidom_loc, nb_subzone;
     /****************************************************
      -----Boucle sous-iteration
     ****************************************************/

        if( nitcfg == 1)
                  {
                     for (E_Int nd = 0; nd < nidom; nd++) //mise a jour metric et vent ale zone cart et 3dhom(3dfull et 2d a la volee)
                         {
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

     // init verrou omp
#include  "FastS/Compute/verrou_lhs_init.cpp"
        //---------------------------------------------------------------------
        // -----Boucle sur num.les domaines de la configuration
        // ---------------------------------------------------------------------
        E_Int shift_zone=0; E_Int shift_wig=0; E_Int shift_coe=0; E_Int nd_current=0;E_Int shift_rk4=0;  E_Int shift_grad=0;
        E_Int shift_a1=0;   E_Int shift_a3=0;  E_Int shift_a4=0;
	E_Float rhs_end=0; E_Int ncells=0;

        E_Int nbtask = ipt_omp[nitcfg-1];
        E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

        for (E_Int ntask = 0; ntask < nbtask; ntask++)
          {
             E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
             E_Int nd = ipt_omp[ pttask ];

             E_Int cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];

             shift_zone=0; shift_wig=0; shift_coe=0; shift_grad=0;
             for (E_Int n = 0; n < nd; n++)
               {
                 shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
                 shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
                 if(param_int[n][ KFLUDOM ]==2){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*3;}
                 if(param_int[n][ KFLUDOM ]==8){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*4;}

                 if(param_int[n][ IFLOW ]==4)
                  {  shift_a1   = shift_a1   + param_int[n][ NDIMDX ]*9;
                     shift_a3   = shift_a3   + param_int[n][ NDIMDX ]*6;
                     shift_a4   = shift_a4   + param_int[n][ NDIMDX ]*param_int[n][ NEQ_LBM ];
                  }

                 if(param_int[n][ITYPZONE ]==4){  shift_grad = shift_grad + param_int[n][ NDIMDX ]*param_int[n][ NEQ ]*3;}
               }

             E_Int lmin = 10;
             if (param_int[nd][ITYPCP] == 2) lmin = 4;

/*
             if (param_int[nd][ITYPZONE] == 4)
              {
//#include       "../../FastP/FastP/Compute/rhs.cpp"
//               shift_grad = shift_grad + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ]*3;
              }
             else{
*/
             if (param_int[nd][IFLOW] == 4 && nitcfg == 1)
              {
                lmin = 4;
                //Quand on fait du couplage : LBM avance en une iteration
//#include        "FastASLBM/Compute/rhs.cpp"
              }
              else if (param_int[nd][IFLOW] != 4)
              {
#include        "FastS/Compute/rhs.cpp"
              }

          } //Fin boucle sur task/zones pour calcul RHS


#ifdef Conservatif
#include   "FastS/Compute/cp_debitIBM.cpp"
#endif

#include   "FastS/Compute/lhs.cpp"
  } // Fin zone // omp


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
  if( omp_mode == 1) { Nbre_thread_actif_loc = 1; ithread_loc = 1;}
  else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

  E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
  E_Int socket          = (ithread-1)/thread_parsock +1;
  E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

  E_Float rhs_end=0;

  E_Int nitcfg_loc= 1;

  E_Int nbtask = ipt_omp[nitcfg_loc-1];
  E_Int ptiter = ipt_omp[nssiter+ nitcfg_loc-1];

  //calcul du sous domaine a traiter par le thread
  for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
        E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
        E_Int nd = ipt_omp[ pttask ];

//#include "FastASLBM/INTERP/stck.cpp"
     }
  } // Fin zone // omp
/*
*/



for (E_Int nd = 0; nd < nidom; nd++)
{   
    if ( param_int[nd][LEVEL] == 1 && nitcfg == 3)
    {
        E_Int pt_interp = param_int[nd][PT_INTERP];
        E_Int nrac = param_int[nd][pt_interp];
        for (E_Int rac=0; rac < nrac; rac++) // Boucle sur les differents raccords
        {
            E_Int pt_racInt = param_int[nd][ pt_interp + rac +1];
            E_Int pos_p1 = param_int[nd][pt_racInt +1];
            pos_p1+=1; if(pos_p1==3){pos_p1=0;}
            param_int[nd][pt_racInt+1] = pos_p1;
        }
    }
}


#ifdef TimeShow
#ifdef _OPENMP
  time_COM = omp_get_wtime();
#endif
#ifdef _OPENMP
  time_init = omp_get_wtime();
#endif
#endif

/*
*/

  //
  //
  //FillGhostcell si mise a jour necessaire et transfer dans C layer
  //
  //
E_Int autorisation_bc[nidom];
for (E_Int nd = 0; nd < nidom; nd++) { autorisation_bc[nd]=1;}

//calcul Sij sur zone NS pour couplage NS_LBM
if(flag_NSLBM==1)
  {
    E_Int nitcfg_loc= 1;
    for (E_Int ip2p = 1; ip2p < param_int_tc[1]+1; ++ip2p)
       {
         E_Int sizecomIBC = param_int_tc[2];
         E_Int sizecomID  = param_int_tc[3+sizecomIBC];
         E_Int shift_graph = sizecomIBC + sizecomID + 3;

         E_Int ech  = param_int_tc[ip2p+shift_graph];
         E_Int dest = param_int_tc[ech];
         if (dest == rank)  // Intra Process
         {
           E_Int numpassage1    = 1;
           E_Int rk1            = param_int[0][RK];
           E_Int exploc1        = param_int[0][EXPLOC];
           E_Int bidim          =1; //on annule les d/dz
           E_Int TypeTransfert1 =0; // uniquement ID

           K_FAST::compute_sij(iptro_CL, iptS, iptvol, param_int_tc, param_real_tc , param_int, param_real, ipt_omp,
                               TypeTransfert1, it_target, nidom, ip2p, bidim, nitcfg_loc, nssiter, rk1, exploc1, numpassage1);
         }
       } //loop ip2p
  }
/*
*/

E_Int nitcfg_stk = nitcfg;
if (lexit_lu ==0 && layer_mode==1)
{
  //remplissage ghost transfert
  #include "FastS/Compute/transfert_multiblock.cpp"
  //E_Int cycl;

  for (E_Int nd = 0; nd < nidom; nd++)
  {
    // flag pour transfert optimiser explicit local
    autorisation_bc[nd]=0;
    if (param_int[0][EXPLOC] == 1)    //if (rk == 3 && exploc == 2)
    {
        cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];

        // modif de nitcfg pour appliquer les BC partout
        if      (nitcfg_stk%cycl == cycl/2 -1       && cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}
        else if (nitcfg_stk%cycl == cycl/2 + cycl/4 && cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}
        else if (nitcfg_stk%cycl == cycl-1          && cycl != 4 ){ nitcfg = 1; autorisation_bc[nd] = 1;}
        else if((nitcfg_stk%cycl == 1 || nitcfg_stk%cycl == cycl/2 || nitcfg_stk%cycl== cycl-1) && cycl == 4 ) { nitcfg = 1; autorisation_bc[nd] = 1; }
    }
    else {autorisation_bc[nd] = 1;}
  }//loop zone
} //lexit



//calcul sequentiel flux IBM et data CLs (lund, wallmodel) a toutes les ssiter
for (E_Int nd = 0; nd < nidom; nd++)
  {
     if( param_int[nd][IFLOW] !=4)
      {
        #include  "FastS/Compute/data_ibm_bc.cpp"
      }
  }

if(lexit_lu ==0 && layer_mode==1)
{
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

   //E_Int Nbre_socket   = NBR_SOCKET;
   E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;
   if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
   else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}


    //remplissage ghost transfert
    #include "Fast/bcs.cpp"
 }//fin zone omp BCs
} //test lexit


  delete [] ipt_timecount;
  delete [] iptro_CL;
  delete [] iptro_ssiter;
  delete [] iptrom_CL;
  delete [] iptroS_CL;

  return ibord_ale;
}
