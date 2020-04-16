# include "fastS.h"
# include "fastP.h"
# include "fastLBM.h"
# include "fastc.h"
# include "fast.h"
# include "../Fast/param_solver.h"
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

E_Int K_FAST::gsdr3( 
  E_Int**& param_int  , E_Float**& param_real ,
  E_Int& nidom        , E_Int& nitrun         , E_Int&  nitcfg    , E_Int&  nitcfg_last, E_Int&  nssiter, E_Int& it_target , E_Int&  first_it,
  E_Int& kimpli       , E_Int& lssiter_verif  , E_Int& lexit_lu   , E_Int& omp_mode    , E_Int& layer_mode, E_Int& mpi,
  E_Int& nisdom_lu_max, E_Int& mx_nidom       , E_Int& ndimt_flt  , E_Int& ndim_grad   ,
  E_Int& threadmax_sdm, E_Int& mx_synchro,
  E_Int& nb_pulse     , 
  E_Float& temps,
  E_Int* ipt_ijkv_sdm  ,
  E_Int* ipt_ind_dm_omp, E_Int* ipt_topology    , E_Int* ipt_ind_CL    , E_Int* ipt_lok,  E_Int* verrou_lhs,  E_Int& vartype, E_Float* timer_omp,
  E_Int*     iptludic  , E_Int*   iptlumax      ,
  E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom  , E_Int** ipt_ng_pe      ,
  E_Float* ipt_VectG      , E_Float* ipt_VectY   , E_Float** iptssor       , E_Float** iptssortmp    , E_Int* ipt_ssor_size, E_Float* ipt_drodmd,
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
  E_Float*   iptroflt        , E_Float*   iptroflt2       , E_Float*  iptgrad      , E_Float*  iptwig       , E_Float*   iptstat_wig  ,
  E_Float*   iptdrodm        , E_Float*   iptcoe          , E_Float*  iptrot       , E_Float**& iptdelta     , E_Float**& iptro_res, E_Float**& iptdrodm_transfer,
  E_Int*&    param_int_tc    , E_Float*& param_real_tc    , E_Int*& linelets_int   , E_Float*& linelets_real , 
  E_Int&     taille_tabs     , E_Float*& stock            , E_Float*& drodmstock   , E_Float*& constk        , E_Float** iptsrc,
  E_Float*   feq             , E_Float**& feq_transfer )


 {
   
#ifdef TimeShow  
  E_Int rank =0; 
  E_Float* ipt_timecount = new E_Float[5];
  ipt_timecount[0:4] = 0.0;

  E_Int nbpointsTot =0;

  for (E_Int nd = 0; nd < nidom; nd++)
      {
        nbpointsTot = nbpointsTot + param_int[nd][ IJKV    ]*param_int[nd][ IJKV +1 ]*param_int[nd][ IJKV +2 ];
      }
#else
  E_Float* ipt_timecount = NULL;        
#endif

#ifdef TimeShow

#ifdef _MPI
  if(mpi)
  {
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);  
  }
#endif

  //ofstream outputfile;
  //std::ostringstream tmp;

  //tmp << "Output" << std::setw(4) << std::setfill('0') << std::to_string(rank);
  //std::string filename = tmp.str();

 // outputfile.open(filename, ios::app);
  time_init = 0;
#ifdef _OPENMP  
  time_init = omp_get_wtime();
#endif
    
#endif


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
      //E_Float** iptro_ssiter;
      E_Float** iptro_ssiter  = new E_Float*[nidom];
      E_Int ishift,lfwmean;
      //E_Float** iptro_CL; 
      //E_Float* iptro_CL[nidom];
      E_Float** iptro_CL  = new E_Float*[nidom];

      E_Int rk =  param_int[0][RK];
      E_Int exploc = param_int[0][EXPLOC];
      E_Int numpassage = 1;

      for (E_Int nd = 0; nd < nidom; nd++)
      {
        if ( param_int[nd][IFLOW]!= 4)
         {
           // J-Scheme + Constantinescu rk2 + rk3 local + rk2local test + Tang & Warnecke
           if( param_int[0][ITYPCP]==2 and ( param_int[0][EXPLOC]== 1 or param_int[0][EXPLOC]== 2 or param_int[0][EXPLOC]== 4 or param_int[0][EXPLOC]== 5) )
	    {
	      if (nitcfg%2 != 0){iptro_ssiter[nd] = iptro[nd];  iptro_CL[nd] = iptrotmp[nd]; ishift  =1; lfwmean  = 1; }
	      else {
	            if (nitcfg != param_int[0][NSSITER]){iptro_ssiter[nd] = iptrotmp[nd]; iptro_CL[nd] = iptro[nd]; ishift  = -1; lfwmean = 0;}
                    else                                {iptro_ssiter[nd] = iptrotmp[nd]; iptro_CL[nd] = iptro[nd];}
                   }

	    }
           else  // Explicit global ou Implicit
	    {
	      if( nitcfg == 1) { iptro_ssiter[nd] = iptro[nd];  iptro_CL[nd] = iptrotmp[nd]; ishift  = 1; lfwmean  = 1; }
	      else
	      {      
               if (param_int[0][ ITYPCP ] < 2                )                                            {iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptrotmp[nd]; ishift =-1; lfwmean  = 0;} // Implicit
               if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==0 && nitcfg != param_int[nd][ NSSITER ] )      {iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptro[nd];    ishift =-1; lfwmean  = 0;} // Explicite
	       else if (param_int[0][ ITYPCP ] == 2 && nitcfg%2==0 && nitcfg == param_int[nd][ NSSITER ] ) {iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptro[nd];                             } // Explicite
               if  (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg != param_int[nd][ NSSITER ])      {iptro_ssiter[nd] = iptro[nd];     iptro_CL[nd] = iptrotmp[nd]; ishift =-1; lfwmean  = 0;} // Explicite
               else if  (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg == param_int[nd][ NSSITER ]) {iptro_ssiter[nd] = iptro[nd];     iptro_CL[nd] = iptrotmp[nd];                          } // Explicite
	      }
	    }
         }//NS ou LBM
        else{ iptro_ssiter[nd] = iptrotmp[nd];  iptro_CL[nd] = iptrotmp[nd]; }

      }//zone

      //
      //Calcul taille tableau ssor par thread et mise a jour Nombre sous_iter pour implicit
      //
      E_Int nd_subzone = 0;
      for (E_Int nd = 0; nd < nidom; nd++)
      {
        //mise a jour Nombre sous_iter pour implicit (simplification gestion OMP)
        if( (param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1 || lssiter_verif ==1 ) && param_int[nd][ ITYPZONE] != 4)
          {
           E_Int* ipt_nisdom_residu   =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter;                //nisdom_residu(nssiter)
           E_Int* ipt_it_bloc         =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter*2;    //it_bloc(nidom)

           if(ipt_nisdom_residu[nitcfg-1] != 0) ipt_it_bloc[0] +=1;
          }

    E_Int* ipt_nidom_loc= ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;


	if( (param_int[ nd ][ NB_RELAX ] > 1 || param_int[ nd ][ LU_MATCH ]==1) && param_int[nd][ ITYPZONE] != 4)
	  {
#ifdef _OPENMP
	    E_Int  Nbre_thread_actif = __NUMTHREADS__;
#else
	    E_Int Nbre_thread_actif = 1;
#endif
	    E_Int nfic_ij = param_int[ nd ][ NIJK + 3 ];
	    E_Int nfic_k  = param_int[ nd ][ NIJK + 4 ];
	    E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
	    E_Int  nb_subzone    = ipt_nidom_loc [nitcfg-1];                                           //nbre sous-zone a la sousiter courante

	    for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
	      {
		if (omp_mode == 1)
		  {
		    E_Int       Ptomp = param_int[nd][PT_OMP];
		    E_Int  PtrIterOmp = param_int[nd][Ptomp];   
		    E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

		    for (E_Int i = 0; i < Nbre_thread_actif; i++)
		      {
			if (param_int[nd][PtZoneomp + i] != - 2)
			  {
			    E_Int* ipt_ind_dm_thread = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (param_int[nd][PtZoneomp + i]) * 6;
			    E_Int indice             = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
		            ipt_ssor_size[indice]=(ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] +1 +2*param_int[nd][ NIJK +3 ])*
	                                	  (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] +1 +2*param_int[nd][ NIJK +3 ])*
		                                  (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] +1 +2*param_int[nd][ NIJK +4 ]);
			  }
                         else { E_Int indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
                                ipt_ssor_size[indice]=0;
                              }
		      }
		  }
		else //ompmode = 0
		  {
                    E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
		    E_Int ipt_topology_socket[3];
		    E_Int ipt_ind_dm_thread[6];
		    E_Int lmin = 10;
		    if (param_int[nd][ITYPCP] == 2) lmin = 4;

		    for (E_Int i = 0; i < Nbre_thread_actif; i++)
		      {
			E_Int ithread = i + 1;
			E_Int indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;

			indice_boucle_lu_(nd, ithread, Nbre_thread_actif, lmin,
					  ipt_ind_dm_loc,
					  ipt_topology_socket, ipt_ind_dm_thread);

			if (ithread > ipt_topology_socket[0]*ipt_topology_socket[1]*ipt_topology_socket[2])
		         { ipt_ssor_size[indice]=0;
			  //break;
                         }
                        else
                         {
		          ipt_ssor_size[indice]=(ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] +1 +2*param_int[nd][ NIJK +3 ])*
	                                        (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] +1 +2*param_int[nd][ NIJK +3 ])*
		                                (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] +1 +2*param_int[nd][ NIJK +4 ]);
                         }
		      }
		  }// ompmode
	      }//loop subzone
	  }  // if relax
      } // loop zone

E_Float deb_calcul;

//cout << "nidom= " << nidom << endl;
//cout << "nitcfg= " << nitcfg << endl; 
//if (nitcfg==4){deb_calcul = omp_get_wtime();}
/****************************************************
----- Debut zone // omp
****************************************************/

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
   E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;
   E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;
   E_Int* ipt_ind_dm_omp_thread  = ipt_ind_dm_socket  + 6;

   E_Int ipt_ind_sdm_thread[6];
   E_Int ipt_ind_coe_thread[6];
   E_Int ipt_ind_grad_thread[6];

   E_Int* ipt_nidom_loc, nb_subzone;
     /****************************************************
      -----Boucle sous-iteration
     ****************************************************/
   //cout << "nticfg= " << nitcfg << endl;
        if( nitcfg == 1)
                  {
                     for (E_Int nd = 0; nd < nidom; nd++) //mise a jour metric et vent ale zone cart et 3dhom(3dfull et 2d a la volee)
                         {
                           if(param_int[nd][LALE]==1 && param_int[nd][ITYPZONE]!= 4 ) //maillage indeformable structure
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

                   // calcul metric si maillage deformable structure
#include           "../../FastS/FastS/Metric/cp_metric.cpp"
                   }


	/// calcul de ndim

	E_Int ndim = 0;
	//cout<<"ndim= "<< ndim << endl;
        //---------------------------------------------------------------------
        // -----Boucle sur num.les domaines de la configuration
        // ---------------------------------------------------------------------
        E_Int shift_zone=0; E_Int shift_wig=0; E_Int shift_coe=0; E_Int nd_current=0;E_Int shift_rk4=0;  E_Int shift_grad=0;
	if  (param_int[0][EXPLOC] == 0 and param_int[0][RK] == 4 )//or param_int[0][EXPLOC] == 0 and param_int[0][RK] == 5)
	  {
	    for (E_Int nd = 0; nd < nidom; nd++)
	      {
		shift_rk4 = shift_rk4 + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
	      }
	  }
        E_Float rhs_end=0;

	shift_rk4 = shift_rk4*(nitcfg - 1);
	//cout << "shift_rk4= " << shift_rk4 << endl;
	

        if (param_int[0][IMPLICITSOLVER] == 1 && layer_mode == 1) { ipt_norm_kry[ithread-1]=0.; }
        //calcul du sous domaine a traiter par le thread 
        for (E_Int nd = 0; nd < nidom; nd++)
          {


	    E_Int cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];
	    //if (param_int[0][EXPLOC] == 2 and param_int[0][RK] == 3 and cycl != 4 and nitcfg%cycl==cycl/4)
	    //{
	    //	 ndim = param_int[0][SHIFTLOCAL];
		 //cout << "ndim= "<< ndim << endl;
	    //		}
	    //else
	    //   {   
	    //	 ndim = 0;
	    //  }
   	    if (param_int[0][EXPLOC] == 2 and param_int[0][RK] == 3 and cycl != 4 and nitcfg%cycl==cycl/4)
	      {
		for (E_Int nd = 0; nd < nidom; nd++)
		  {
		    ndim = ndim + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
		  }

	      }
	    else
	       {   
	    	 ndim = 0;
	       }
   
           E_Int lmin = 10;
           if (param_int[nd][ITYPCP] == 2) lmin = 4;



          if (param_int[nd][ITYPZONE] == 4)
           {
             if ( ithread == 1)
             {
#include       "../../FastP/FastP/Compute/rhs.cpp"
               shift_grad = shift_grad + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ]*3;
             }
           }
          else{
              if (param_int[nd][IFLOW] == 4){
               lmin =1;
#include       "../../FastLBM/FastLBM/Compute/rhs.cpp"
              }
              else{
#include       "../../FastS/FastS/Compute/rhs.cpp"
              }
           }
          shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
          shift_wig  = shift_wig  + param_int[nd][ NDIMDX ]*3;
          shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];

	  //E_Float fin_zone = omp_get_wtime();

	  //if(ithread==1){cout <<"zone : "<< nd <<" "<< "temps= " << fin_zone - deb_zone <<" "<<"cycle =  " << cycl << endl;}
	   

          } //Fin boucle sur zones pour calcul RHS
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
          //if(ithread==1) printf(" time rhs= %g %g  %d %d  \n",timer_omp[ cpu_perthread], rhs_end - rhs_begin , cpu_perthread, nitcfg);

#include   "../../FastS/FastS/Compute/lhs.cpp"
          //
          //finalisation timer pour omp "dynamique"
          //
#ifdef _OPENMP  
          E_Float     lhs_end = omp_get_wtime();
#else  
          E_Float     lhs_end = 0.;
#endif
          timer_omp[ cpu_perthread +1 ] += lhs_end- rhs_end;
          // if(ithread==1) printf(" time lhs= %g \n",lhs_end - rhs_end );
          

} // Fin zone // omp



#ifdef TimeShow
#ifdef _OPENMP  
     time_COM = omp_get_wtime();
#endif
//     outputfile << "Time in compute (gsdr3 omp//) " << time_COM - time_init << std::endl;
//     outputfile << "Time adim: " <<  (time_COM - time_init)/nbpointsTot << " nidom " << nidom << std::endl;
     // std::cout << " sizes (ni,nj,nk):" << std::endl;
     // for (E_Int nd = 0; nd < nidom; nd++)
     //      {
     //      std::cout << param_int[nd][ IJKV    ] << "," << param_int[nd][ IJKV +1 ] << "," << param_int[nd][ IJKV  +2  ] << std::endl;
     //      } 
//     outputfile.close();
#ifdef _OPENMP  
     time_init = omp_get_wtime();
#endif 
#endif 


 
  //
  //
  //FillGhostcell si mise a jour necessaire et transfer dans C layer 
  //
  //
if(lexit_lu ==0 && layer_mode==1)
{   

  //Swap (call to setInterpTransfer)

  //E_Float trans_begin = omp_get_wtime();


E_Int cycl;
E_Float deb_dtlocal;
E_Float tmps;
    
  if (param_int[0][EXPLOC] == 0)
    {
      K_FASTC::setInterpTransfersFast(iptro_CL, vartype, param_int_tc, param_real_tc , param_int, param_real,
                                   linelets_int, linelets_real, it_target, nidom, ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage);
    }//test exploc==0

 
  if (param_int[0][EXPLOC] == 2 and param_int[0][RK] == 3) 
   {
     K_FASTS::dtlocal2para_c(iptro, iptrotmp, param_int_tc, param_real_tc, param_int, param_real, iptdrodm, iptcoe, stock, drodmstock, constk, nitcfg, omp_mode, taille_tabs, nidom);

     for (E_Int nd=0; nd < nidom; nd++)
       {
	 cycl = nssiter/param_int[nd][LEVEL];
	 if (nitcfg%cycl == 0 and nitcfg != nssiter)
	   {
	     E_Float* ptsave  = iptro[nd]; 
	     iptro[nd] = iptrotmp[nd]; 
	     iptrotmp[nd] = ptsave;
	   } 
       }

      K_FASTC::setInterpTransfersFast(iptro_CL, vartype, param_int_tc, param_real_tc , param_int, param_real,
                             linelets_int, linelets_real, it_target   , nidom        , ipt_timecount, mpi       , nitcfg, nssiter, rk, exploc, numpassage);


// achtung
//     recup3para_c(iptro, param_int_tc, param_real_tc, param_int, stock, nitcfg, omp_mode, taille_tabs, nidom);

     numpassage=2;
     if (nitcfg%2==0)
      {
        K_FASTC::setInterpTransfersFast(iptro_CL, vartype, param_int_tc, param_real_tc , param_int, param_real,
                               linelets_int, linelets_real, it_target   , nidom        , ipt_timecount, mpi       , nitcfg, nssiter, rk, exploc, numpassage);
      }
 

   } // fin boucle test dtlocal

  


  //
  //E_Float trans_end = omp_get_wtime();
  //E_Float trans_duree = trans_end - trans_begin;
  //cout << "tps_trans  : "<< trans_duree  << endl;
  //
  //
  //Apply BC (parcour Zones) + reinitialisation verrou pour calcul rhs
  //
  //
    
 E_Int autorisation_bc[nidom];
 E_Int nitcfg_stk = nitcfg;
 
for (E_Int nd = 0; nd < nidom; nd++)
  {
   // flag pour transfert optimiser explicit local
   autorisation_bc[nd]=0;
   if (rk == 3 and exploc == 2) 
      {
	cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];
	  
	if      (nitcfg_stk%cycl == cycl/2 -1       and cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}
	else if (nitcfg_stk%cycl == cycl/2 + cycl/4 and cycl != 4) { nitcfg = 1; autorisation_bc[nd] = 1;}	    
	else if (nitcfg_stk%cycl == cycl-1          and cycl != 4 ){ nitcfg = 1; autorisation_bc[nd] = 1;}
	else if((nitcfg_stk%cycl == 1 or nitcfg_stk%cycl == cycl/2  or nitcfg_stk%cycl== cycl-1) and cycl == 4 ) { nitcfg = 1; autorisation_bc[nd] = 1; }
      } 
    else {autorisation_bc[nd] = 1;}

    if(nitcfg==nitcfg_last-1 and param_int[nd][ITYPZONE] !=4  and param_int[nd][IFLOW] !=4)
    {
     ///mise a jour moyenne plan Lund si necessaire
     E_Int pt_bcs = param_int[nd][PT_BC];
     E_Int nb_bc  = param_int[nd][ pt_bcs ];
     E_Float* ipt_data=NULL;
     for ( E_Int ndf = 0; ndf < nb_bc; ndf++ )
     {
      E_Int pt_bc  = param_int[nd][pt_bcs+ 1 + ndf];

      E_Int idir   = param_int[nd][pt_bc + BC_IDIR];
      E_Int nbdata = param_int[nd][pt_bc + BC_NBDATA];
      E_Int bc_type= param_int[nd][pt_bc + BC_TYPE];

      if(bc_type==19) 
      { 
         E_Int* iptsize_data = param_int[nd] + pt_bc + BC_NBDATA + 1;
         E_Int* ind_fen      = param_int[nd] + pt_bc + BC_FEN;
         E_Int  inc_bc[3];
         if ( idir <= 2 ) 
         { inc_bc[0] = ind_fen[3] - ind_fen[2] + 1; // nombre element de la fenetre dans la direction J
           inc_bc[1] = ind_fen[2]; // debut indice j
           inc_bc[2] = ind_fen[4]; // debut indice k
         }  
         else if ( idir <= 4 )
         { inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
           ind_fen[0]; // debut indice i
           inc_bc[2] = ind_fen[4]; // debut indice k
         }  
         else
         {inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
          inc_bc[1] = ind_fen[0]; // debut indice i
          inc_bc[2] = ind_fen[2]; // debut indice j
         }

         if ( nbdata != 0 ) ipt_data = param_real[nd] + param_int[nd][pt_bcs + 1 + ndf + nb_bc];

         E_Float* iptAvgPlanLund = ipt_data + 5*iptsize_data[0];
         E_Float* iptParamLund   = ipt_data + 10*iptsize_data[0];

         E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

         mj_lund_planrecyl_(nd, idir, param_int[nd] , ind_fen, ipt_ind_dm_loc, inc_bc, iptsize_data[0],
                            iptParamLund,  iptro[nd], iptAvgPlanLund); 
      }//lund
     }//ndf
    }//nit_last
  }//loop zone
 
  E_Int ndim=0;
  //return ndim;

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


   E_Int* ipt_ind_dm_thread; E_Int* ipt_topology_socket; E_Int* ipt_ind_dm_socket;
   E_Int nd_current =0;
   for (E_Int nd = 0; nd < nidom; nd++)
     {
       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
       E_Int  nb_subzone    = ipt_nidom_loc [nitcfg-1];                                           //nbre sous-zone a la sousiter courante

       E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
       E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
       E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
       E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;

       for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
	 {
	   E_Int ndo   = nd;

	   E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

	   if (omp_mode == 1)
	     { 
	       E_Int       Ptomp = param_int[nd][PT_OMP];
	       E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
	       E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

	       Nbre_thread_actif_loc = param_int[nd][ PtZoneomp  + Nbre_thread_actif ];
	       ithread_loc           = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
	       ipt_ind_dm_thread     = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;

	       if (ithread_loc == -1) { continue;}
	     }
	   else
	     { 
	       ipt_topology_socket = ipt_topology       + (ithread-1)*3;
	       ipt_ind_dm_socket   = ipt_ind_dm_omp     + (ithread-1)*12;
	       ipt_ind_dm_thread   = ipt_ind_dm_socket  +6;

	       E_Int lmin = 10;
	       if (param_int[nd][ITYPCP] == 2) lmin = 4;

               indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, lmin,
                                 ipt_ind_dm_loc,
                                 ipt_topology_socket, ipt_ind_dm_thread);
	     }

	   if (autorisation_bc[nd] == 1)
	    {
               if (param_int[nd][ITYPZONE ] !=4) // zone structuree
                {
                  if(param_int[nd][IFLOW] != 4)
                    {
	              E_Int ierr = K_FASTS::BCzone(nd, lrhs , nitcfg_stk, lcorner, param_int[nd], param_real[nd], npass,
                                                   ipt_ind_dm_loc, ipt_ind_dm_thread, 
                                                   ipt_ind_CL_thread, ipt_ind_CL119,  ipt_ind_CLgmres, ipt_shift_lu,
                                                   iptro_CL[nd] , ipti[nd]            , iptj[nd]    , iptk[nd]       ,
                                                   iptx[nd]     , ipty[nd]            , iptz[nd]    ,
				                   iptventi[nd] , iptventj[nd]        , iptventk[nd], iptro_CL[nd]);

	              correct_coins_(nd,  param_int[nd], ipt_ind_dm_thread , iptro_CL[nd]);
                    }
                   else
                    {
                     E_Int ierr = K_FASTLBM::BCzone(nd, lrhs , nitcfg_stk, lcorner, param_int[nd], param_real[nd], npass,
                                                    ipt_ind_dm_loc, ipt_ind_dm_thread, 
                                                    ipt_ind_CL_thread, ipt_ind_CL119,  ipt_ind_CLgmres, ipt_shift_lu,
                                                    iptro_CL[nd] , ipti[nd]            , iptj[nd]    , iptk[nd]       ,
                                                    iptx[nd]     , ipty[nd]            , iptz[nd]    ,
                                                    iptventi[nd] , iptventj[nd]        , iptventk[nd], iptro_CL[nd]); 
                    }
                }
               else
                {
                     E_Int lrhs=0; E_Int lcorner=0;
                     if (ithread ==1) K_FASTP::BCzone( nd, lrhs, lcorner,
                                                        param_int[nd], param_real[nd],
                                                        npass, temps,
                                                        ipt_ind_CL119  , 
                                                        ipt_ng_pe[nd]          , iptro_CL[nd]   ,
                                                        ipti[nd] , iptventi[nd], iptx[nd] , ipty[nd] , iptz[nd] );
                }
	    }//autorisation

           //Reinitialisation verrou omp
           //
           E_Int l =  nd_current*mx_synchro*Nbre_thread_actif  + (ithread_loc-1)*mx_synchro;
           for (E_Int i = 0;  i < mx_synchro ; i++) { ipt_lok[ l + i ]  = 0; }
	   nd_current +=1;

	 }//loop souszone
     }//loop zone

 }//fin zone omp


  nitcfg = nitcfg_stk;

//     cout << "nidom= "<< nidom << endl;	 

   //cout << "coucou" << endl;



 }//test exit_lu


    //
    //omp "dynamic" balance
    //
//    if(nitrun%5==0)
//    {  
//#include "HPC_LAYER/OPTIMIZE_DISTRIB_OMP.h"
//
//    }


   delete [] ipt_timecount; delete [] iptro_CL;  delete [] iptro_ssiter; 
   return ibord_ale;

}
