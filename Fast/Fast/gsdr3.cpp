# include "FastS/fastS.h"
//# include "fastP.h"
# include "FastLBM/fastLBM.h"
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

E_Int K_FAST::gsdr3(E_Int**& param_int          , E_Float**& param_real       ,
		    E_Int& nidom                , E_Int& nitrun               ,
		    E_Int&  nitcfg              , E_Int&  nitcfg_last         ,
		    E_Int&  nssiter             , E_Int& it_target            ,
		    E_Int&  first_it            , E_Int& kimpli               ,
		    E_Int& lssiter_verif        , E_Int& lexit_lu             ,
		    E_Int& omp_mode             , E_Int& layer_mode           ,
		    E_Int& mpi                  , E_Int& nisdom_lu_max        ,
		    E_Int& mx_nidom             , E_Int& ndimt_flt            ,
		    E_Int& threadmax_sdm        , E_Int& mx_synchro           ,
		    E_Int& nb_pulse             , E_Float& temps              ,
		    E_Int* ipt_ijkv_sdm         , E_Int* ipt_ind_dm_omp       ,  E_Int* ipt_omp           ,
		    E_Int* ipt_topology         , E_Int* ipt_ind_CL           ,
		    E_Int* ipt_lok              , E_Int* verrou_lhs           ,
		    E_Int& vartype              , E_Float* timer_omp          ,
		    E_Int*     iptludic         , E_Int*   iptlumax           ,
		    E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom     ,
		    E_Int** ipt_ng_pe           , E_Int** ipt_nfconn          ,
		    E_Int** ipt_nfindex         , E_Float* ipt_VectG          ,
		    E_Float* ipt_VectY          , E_Float** iptssor           ,
		    E_Float** iptssortmp        , E_Int* ipt_ssor_size        ,
		    E_Float* ipt_drodmd         , E_Float* ipt_Hessenberg     ,
		    E_Float** iptkrylov         , E_Float** iptkrylov_transfer,
		    E_Float* ipt_norm_kry       , E_Float** ipt_gmrestmp      ,
		    E_Float* ipt_givens         , E_Float*   ipt_cfl          ,
		    E_Float**  iptx             , E_Float**  ipty             ,
		    E_Float**    iptz           , E_Float**  iptCellN         ,
		    E_Float**  iptCellN_IBC     , E_Int** ipt_degen           ,
		    E_Float**& iptro            , E_Float**& iptro_m1         ,
		    E_Float**&  iptrotmp        , E_Float**& iptrof           ,
		    E_Float**  iptmut           , E_Float*   ipt_mutd         ,
		    E_Float**  ipti             , E_Float**  iptj             ,
		    E_Float** iptk              , E_Float** iptvol            ,
		    E_Float**  ipti0            , E_Float**  iptj0            ,
		    E_Float** iptk0             , E_Float**  ipti_df          ,
		    E_Float**  iptj_df          , E_Float** iptk_df           ,
		    E_Float** iptvol_df         , E_Float**  iptventi         ,
		    E_Float**  iptventj         , E_Float** iptventk          ,
		    E_Float**& iptrdm           , E_Float*   iptroflt         ,
		    E_Float*   iptroflt2        , E_Float*  iptgrad           ,
		    E_Float*  iptwig            , E_Float*   iptstat_wig      ,
		    E_Float* iptflu             , E_Float*   iptdrodm         ,
		    E_Float*   iptcoe           , E_Float*  iptrot            ,
		    E_Float**& iptdelta         , E_Float**& iptro_res        ,
		    E_Float**& iptdrodm_transfer, E_Int*&    param_int_tc     ,
		    E_Float*& param_real_tc     , E_Int*& linelets_int        ,
		    E_Float*& linelets_real     , E_Int&     taille_tabs      ,
		    E_Float*& stock             , E_Float*& drodmstock        ,
		    E_Float*& constk            , E_Float** iptsrc            ,
		    E_Float*   feq              , E_Float**& feq_transfer     ,
		    E_Float**&  iptSpongeCoef   , E_Float**& iptmacro_m1      ,
		    E_Float**& iptdist2ibc      , E_Float**& iptQeq           ,
		    E_Float**&  iptQneq         , E_Float**& iptcellN_IBC_LBM  )
{

#ifdef TimeShow  
  E_Int rank =0; 
  E_Float* ipt_timecount = new E_Float[5];
  ipt_timecount[0:4] = 0.0;

  E_Int nbpointsTot =0;

  for (E_Int nd = 0; nd < nidom; nd++){
      nbpointsTot = nbpointsTot + param_int[nd][ IJKV    ]*param_int[nd][ IJKV +1 ]*param_int[nd][ IJKV +2 ];
  }
#else
  E_Float* ipt_timecount = NULL;        
#endif

#ifdef TimeShow

#ifdef _MPI
  if(mpi){
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

  E_Float** iptmacro=iptro;
  E_Float** iptQ    =iptrotmp; 
  E_Float** iptQstar=iptro_m1;
  
  E_Int ishift,lfwmean;
  E_Float** iptro_ssiter   = new E_Float*[nidom];
  E_Float** iptro_CL       = new E_Float*[nidom];
  E_Float** iptro_CL_macro = new E_Float*[nidom];
  E_Float** iptro_CL_Qneq  = new E_Float*[nidom];
  E_Float** iptro_CL_Qstar = new E_Float*[nidom];
      
  E_Int rk     =  param_int[0][RK];
  E_Int exploc = param_int[0][EXPLOC];
  E_Int numpassage = 1;

  for (E_Int nd = 0; nd < nidom; nd++){
    //NS ou LBM
    if ( param_int[nd][IFLOW]!= 4){
      // explicit local
      if(param_int[0][EXPLOC]== 1 and param_int[0][ITYPCP]==2){
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
	  if  (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg != param_int[nd][ NSSITER ])      {iptro_ssiter[nd] = iptro[nd];     iptro_CL[nd] = iptrotmp[nd]; ishift =-1; lfwmean  = 0;} // Explicite
	  else if  (param_int[0][ ITYPCP ] == 2 && nitcfg%2==1 && nitcfg == param_int[nd][ NSSITER ]) {iptro_ssiter[nd] = iptro[nd];     iptro_CL[nd] = iptrotmp[nd];                          } // Explicite
	}
      }
    }
    else{ // param_int[nd][IFLOW]== 4
      iptro_ssiter[nd]   = iptrotmp[nd];
      iptro_CL[nd]       = iptrotmp[nd];
      iptro_CL_macro[nd] = iptro[nd];
      iptro_CL_Qneq[nd]  = iptQeq[nd];
      iptro_CL_Qstar[nd] = iptro_m1 [nd];
    }
  }//zone

  //
  //Calcul taille tableau ssor par thread et mise a jour Nombre sous_iter pour implicit
  //
#include "prep_NS_implicit.cpp"

  //E_Float deb_calcul;
  //cout << "nidom= " << nidom << endl;
  //cout << "nitcfg= " << nitcfg << endl; 
  //if (nitcfg==4){deb_calcul = omp_get_wtime();}
  
  /****************************************************
----- Debut zone // omp
  ****************************************************/
#include "compute_omp_zone.cpp"

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
# include "fast_c_layer.cpp" 

  //
  //omp "dynamic" balance
  //
  //    if(nitrun%5==0)
  //    {  
  //#include "HPC_LAYER/OPTIMIZE_DISTRIB_OMP.h"
  //
  //    }


  delete [] ipt_timecount;
  delete [] iptro_CL;
  delete [] iptro_ssiter;
  delete [] iptro_CL_macro;
  delete [] iptro_CL_Qneq;
  delete [] iptro_CL_Qstar;
  
  return ibord_ale;
}
