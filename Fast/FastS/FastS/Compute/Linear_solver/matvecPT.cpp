/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
// Mettre a 1 pour un CPU timer
#define TIMER 0

# include "FastC/fastc.h"
# include "FastS/fastS.h"
# include "FastS/param_solver.h"
# include <string.h>
//# include <omp.h>
using namespace std;
using namespace K_FLD;
#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::_matvecPT(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* work; 

  E_Int nitrun, no_vect_test; 
  
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOllO" , &zones , &metrics, &nitrun, &no_vect_test,  &work)) return NULL;
#else 
  if (!PyArg_ParseTuple(args, "OOiiO" , &zones , &metrics, &nitrun, &no_vect_test,  &work)) return NULL;
#endif

  //* tableau pour stocker dimension sous-domaine omp *//
  E_Int threadmax_sdm  = __NUMTHREADS__;

  PyObject* tmp = PyDict_GetItemString(work,"MX_SSZONE");      E_Int mx_sszone       = PyLong_AsLong(tmp);  
            tmp = PyDict_GetItemString(work,"MX_SYNCHRO");     E_Int mx_synchro      = PyLong_AsLong(tmp);  
            tmp = PyDict_GetItemString(work,"MX_OMP_SIZE_INT");E_Int mx_omp_size_int  = PyLong_AsLong(tmp);  
            tmp = PyDict_GetItemString(work,"FIRST_IT");       E_Int first_it        = PyLong_AsLong(tmp);
            tmp = PyDict_GetItemString(work,"mpi");            E_Int mpi             = PyLong_AsLong(tmp);


  PyObject* dtlocArray  = PyDict_GetItemString(work,"dtloc"); FldArrayI* dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc  = dtloc->begin();
  E_Int nssiter = iptdtloc[0];
  E_Int omp_mode = iptdtloc[8];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;


  E_Int lcfl= 0;
 // 
 // 
 // 
 // partie code specifique au Transfer Data
 // 
 // 
 //
 PyObject* pyParam_int_tc;PyObject* pyParam_real_tc; PyObject* iskipArray;
 PyObject* pyLinlets_int; PyObject* pyLinlets_real; 
 FldArrayI* param_int_tc;  E_Int*  ipt_param_int_tc ; FldArrayI* iskip_lu;
 FldArrayF* param_real_tc;
 FldArrayF* linelets_real; FldArrayI* linelets_int;
 E_Int* ipt_iskip_lu;
 E_Float* ipt_param_real_tc;
 E_Float* ipt_linelets_real; E_Int* ipt_linelets_int;

 E_Int lssiter_loc; E_Int lexit_lu;  E_Int lssiter_verif;
 E_Int it_target = iptdtloc[4];

 pyParam_int_tc = PyDict_GetItemString(work,"param_int_tc");
 pyParam_real_tc= PyDict_GetItemString(work,"param_real_tc"); 

 if (pyParam_int_tc != Py_None)
 { K_NUMPY::getFromNumpyArray(pyParam_int_tc , param_int_tc , true); ipt_param_int_tc = param_int_tc -> begin(); }
 else{ ipt_param_int_tc = NULL; }

 if (pyParam_real_tc != Py_None)
 { K_NUMPY::getFromNumpyArray(pyParam_real_tc, param_real_tc, true); ipt_param_real_tc= param_real_tc-> begin(); }
 else{ ipt_param_real_tc = NULL; }

 pyLinlets_int = PyDict_GetItemString(work,"linelets_int");
 if (pyLinlets_int != Py_None)
 {K_NUMPY::getFromNumpyArray(pyLinlets_int, linelets_int, true); ipt_linelets_int = linelets_int->begin();}
 else{ipt_linelets_int = NULL;}  

 pyLinlets_real = PyDict_GetItemString(work,"linelets_real");
 if (pyLinlets_real != Py_None)
 {K_NUMPY::getFromNumpyArray(pyLinlets_real, linelets_real, true); ipt_linelets_real = linelets_real->begin();}
 else{ipt_linelets_real = NULL;}

// Fin partie code specifique au Transfer Data

  E_Int nidom    = PyList_Size(zones);

  FldArrayI n0_flt(nidom); E_Int* ipt_n0_flt = n0_flt.begin();


  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; E_Int** ipt_it_lu_ssdom;
  
  E_Float** ipt_param_real  ;
  E_Float** iptx;          E_Float** ipty;         E_Float** iptz;    E_Float** iptro; E_Float** iptro_m1; E_Float** iptro_p1; E_Float** iptmut;
  E_Float** ipti;          E_Float** iptj;         E_Float** iptk;    E_Float** iptvol;
  E_Float** ipti0;         E_Float** iptj0;        E_Float** iptk0;
  E_Float** ipti_df;       E_Float** iptj_df;      E_Float** iptk_df ; E_Float** iptvol_df;
  E_Float** iptdist;       E_Float** iptventi;     E_Float** iptventj; E_Float** iptventk;
  E_Float** iptrdm;        E_Float** iptkrylov;    E_Float** iptkrylov_transfer;
  E_Float** iptro_sfd;     E_Float** iptdelta;     E_Float** iptfd; E_Float** iptro_zgris; E_Float** iptro_res;
  E_Float** iptCellN  ;    E_Float** iptCellN_IBC; E_Float** ipt_cfl_zones; 
  E_Float** iptssor;       E_Float** iptssortmp;
  E_Float** ipt_gmrestmp;  E_Float** iptdrodm_transfer;

  ipt_param_int     = new E_Int*[nidom*3];
  ipt_ind_dm        = ipt_param_int   + nidom;
  ipt_it_lu_ssdom   = ipt_ind_dm      + nidom;

  iptx              = new E_Float*[nidom*35];
  ipty              = iptx               + nidom;
  iptz              = ipty               + nidom;
  iptro             = iptz               + nidom;
  iptro_m1          = iptro              + nidom;
  iptro_p1          = iptro_m1           + nidom;
  iptro_sfd         = iptro_p1           + nidom;
  iptdelta          = iptro_sfd          + nidom;
  iptro_res         = iptdelta           + nidom;
  iptmut            = iptro_res          + nidom;
  ipti              = iptmut             + nidom;
  iptj              = ipti               + nidom;
  iptk              = iptj               + nidom;
  iptvol            = iptk               + nidom;
  ipti0             = iptvol             + nidom;
  iptj0             = ipti0              + nidom;
  iptk0             = iptj0              + nidom;
  ipti_df           = iptk0              + nidom;
  iptj_df           = ipti_df            + nidom;
  iptk_df           = iptj_df            + nidom;
  iptvol_df         = iptk_df            + nidom;
  iptventi          = iptvol_df          + nidom;
  iptventj          = iptventi           + nidom;
  iptventk          = iptventj           + nidom;
  iptrdm            = iptventk           + nidom;
  iptCellN          = iptrdm             + nidom; 
  ipt_param_real    = iptCellN           + nidom;
  iptkrylov         = ipt_param_real     + nidom;
  iptkrylov_transfer= iptkrylov          + nidom;
  ipt_gmrestmp      = iptkrylov_transfer + nidom;
  iptssor           = ipt_gmrestmp       + nidom;   //ndimdx+2rangee ghost par thread
  iptssortmp        = iptssor            + nidom;   //ndimdx
  ipt_cfl_zones     = iptssortmp         + nidom;   //3composants: attention a la prochaine addition
  iptdrodm_transfer = ipt_cfl_zones      + nidom;
  iptCellN_IBC      = iptdrodm_transfer  + nidom;

  vector<PyArrayObject*> hook;
  PyObject* ssorArray = PyDict_GetItemString(work,"ssors");
    
  E_Int nb_pulse      = 0;                  // =1 => pulse acoustic centree a l'origine
  E_Int ndimdx_max    = 0;
  E_Int ndimdx_sa     = 0;
  E_Int nisdom_lu_max = 0;
  E_Int ndimt         = 0;
  E_Int ndim_drodm    = 0;
  E_Int ndimt_flt     = 0;
  E_Int kles          = 0;
  E_Int kimpli        = 0;
  E_Int kfiltering    = 0;
  E_Int neq_max       = 5;
  E_Int iorder_flt    =10;
  E_Float temps;

  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    // check zone
    PyObject* zone = PyList_GetItem(zones, nd); // domaine i

    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);
                       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(o, hook);

    temps = (nitrun-1)*ipt_param_real[nd][ DTC]; 

    PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i

    if( ipt_param_int[nd][ NDIMDX ] > ndimdx_max ) ndimdx_max = ipt_param_int[nd][ NDIMDX ];
    
    //
    //
    //Pointeur maillage
    //
    //
    if(ipt_param_int[nd][ LALE ]== 0){ GET_XYZ( "GridCoordinates"     , zone, iptx[nd], ipty[nd], iptz[nd])}
    else                             { GET_XYZ( "GridCoordinates#Init", zone, iptx[nd], ipty[nd], iptz[nd])}

    //
    //Pointeur var primitive + visco + distance paroi + cellN
    //
    //
    PyObject* sol_center; PyObject* t;
    sol_center   = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
    t            = K_PYTREE::getNodeFromName1(sol_center, "Density");
    iptro[nd]    = K_PYTREE::getValueAF(t, hook);
    t            = K_PYTREE::getNodeFromName1(sol_center, "Density_M1");
    iptro_m1[nd] = K_PYTREE::getValueAF(t, hook);
    t            = K_PYTREE::getNodeFromName1(sol_center, "Density_P1");
    iptro_p1[nd] = K_PYTREE::getValueAF(t, hook);
    t            = K_PYTREE::getNodeFromName1(sol_center, "cellN");
    if (t == NULL) iptCellN[nd] = NULL;
    else iptCellN[nd] = K_PYTREE::getValueAF(t, hook);

    if(ipt_param_int[nd][ IBC ]== 1){t=K_PYTREE::getNodeFromName1(sol_center, "cellN_IBC"); iptCellN_IBC[nd] = K_PYTREE::getValueAF(t, hook); } 
    else { iptCellN_IBC[nd] = NULL;}

    //Pointeur Vect Krylov
    iptkrylov[nd] = NULL;
    if (ipt_param_int[nd][ IMPLICITSOLVER ] == 1) // 1 = gmres
    {
      t            = K_PYTREE::getNodeFromName1(sol_center, "Kry_0_Density");
      iptkrylov[nd]= K_PYTREE::getValueAF(t, hook);
    }

    iptssor[nd] = NULL;
    if (ipt_param_int[nd][ NB_RELAX ] > 1) // 1 = LUSSOR
      iptssor[nd] = K_NUMPY::getNumpyPtrF(PyList_GetItem(ssorArray, nd));
  
    //Pointeur sfd
    if (ipt_param_int[nd][ SFD ] == 1)
    {
      t            = K_PYTREE::getNodeFromName1(sol_center, "Density_f");
      iptro_sfd[nd]= K_PYTREE::getValueAF(t, hook);
    }
    //Pointeur debug DES: delta et fd compact
    if (ipt_param_int[nd][ SA_INT+ SA_IDES-1 ] >= 2 && ipt_param_int[nd][ SA_DEBUG ]==1)
    {
      t            = K_PYTREE::getNodeFromName1(sol_center, "delta");
      iptdelta[nd] = K_PYTREE::getValueAF(t, hook);
    }    
    //Pointeur extraction RHS
    if (ipt_param_int[nd][ EXTRACT_RES ] == 1)
    {
      t            = K_PYTREE::getNodeFromName1(sol_center, "Res_Density");
      iptro_res[nd]= K_PYTREE::getValueAF(t, hook);
    }

    //Pointeur visqeux: mut, dist, zgris sont en acces compact
    if(ipt_param_int[nd][ IFLOW ] > 1)
      { t          = K_PYTREE::getNodeFromName1(sol_center, "ViscosityEddy");
        iptmut[nd] = K_PYTREE::getValueAF(t, hook);

        if (ipt_param_int[nd][ IFLOW ] ==3) 
            { neq_max = 6;
              //Recherche du domaine RANS/DES le + gros
              if (ipt_param_int[nd][ ILES ] == 0) { if( ipt_param_int[nd][ NDIMDX ] > ndimdx_sa ) ndimdx_sa = ipt_param_int[nd][ NDIMDX ]; }
            }
      }
    else {iptmut[nd] = iptro[nd];}


    //
    //
    // Pointeur metric
    //
    E_Float* dummy;
    GET_TI(METRIC_TI  , metric, ipt_param_int[nd], ipti [nd]  , iptj [nd]  , iptk [nd]  , iptvol[nd]   )
    GET_TI(METRIC_TI0 , metric, ipt_param_int[nd], ipti0[nd]  , iptj0[nd]  , iptk0[nd]  , dummy )
    GET_TI(METRIC_TIDF, metric, ipt_param_int[nd], ipti_df[nd], iptj_df[nd], iptk_df[nd], iptvol_df[nd])

    GET_VENT( METRIC_VENT, metric, ipt_param_int[nd], iptventi[nd], iptventj[nd], iptventk[nd] )

    iptrdm[ nd ]          =  K_NUMPY::getNumpyPtrF( PyList_GetItem(metric, METRIC_RDM ) );
    ipt_ind_dm[ nd ]      =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );
    ipt_it_lu_ssdom[ nd ] =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_ITLU) );

    if( ipt_param_int[ nd ][ MXSSDOM_LU ] > nisdom_lu_max) nisdom_lu_max = ipt_param_int[ nd ][ MXSSDOM_LU ];

    ndimt      += ipt_param_int[nd][ NDIMDX ];
    ndim_drodm += ipt_param_int[nd][ NDIMDX ]*ipt_param_int[nd][ NEQ ];

    if (ipt_param_int[nd][ ITYPCP ]      <= 1     ) kimpli     = 1;
    if (ipt_param_int[nd][ IFLAGFLT ]    == 1     ) kfiltering = 1;
    if (ipt_param_int[nd][ IFLOW ] >= 2 &&  ipt_param_int[nd][ ILES ] == 1 && ipt_param_int[nd][ IJKV+ 2] > 1) kles = 1;

     //Recherche domaine a filtrer
     if (ipt_param_int[nd][ IFLAGFLT ] == 1) 
          { kfiltering     = 1; 
            ipt_n0_flt[nd] = ndimt_flt;
            if (ipt_param_int[nd][ ITYPZONE] == 3)ndimt_flt=ndimt_flt+ (ipt_param_int[nd][ IJKV ]+iorder_flt)*(ipt_param_int[nd][ IJKV +1]+iorder_flt); 
            else     ndimt_flt=ndimt_flt+ (ipt_param_int[nd][ IJKV ]+iorder_flt)*(ipt_param_int[nd][ IJKV +1]+iorder_flt)*(ipt_param_int[nd][ IJKV +2]+iorder_flt);
          }
  } // boucle zone
  
//  // Reservation tableau travail temporaire pour calcul du champ N+1

  E_Int vartype = 2;
  if(neq_max == 6){ vartype = 21;}


  /// Tableau pour filtrage Visbal
  if (kfiltering == 0) neq_max = 0;
  FldArrayF  roflt(ndimt_flt*neq_max) ; E_Float* iptroflt  = roflt.begin();
  FldArrayF  roflt2(ndimt_flt*neq_max); E_Float* iptroflt2 = roflt2.begin();    //optimisation memoire possible pour roflt2 (voir Funk)

  /// Tableau pour  pour modele LES
  ///sauvegarde pointeur pour allouer memoire SGS qu'on peut ecraser une fois calculer mut
  E_Int neq_les  = 0;
  if(kles==1) neq_les = 3;
  FldArrayF  rot(ndimt*neq_les);
  E_Float* iptrot = rot.begin();
  E_Int kwig_stat    = 0;  
  E_Int neq_wig_stat = 0;
  if (kwig_stat ==  1) neq_wig_stat = 3; 
  FldArrayF  stat_wig(ndimt*neq_wig_stat); E_Float* iptstat_wig = stat_wig.begin();

  // INIT du LUSSORTMP
  E_Int size_tot = 0;
  E_Int* ipt_ssor_size = NULL;

  if (ipt_param_int[0][NB_RELAX] > 1 || (ipt_param_int[0][IMPLICITSOLVER] == 1 && ipt_param_int[0][NB_RELAX] != 0))
    for (E_Int nd = 0; nd < nidom; nd++)
      size_tot += ipt_param_int[nd][NDIMDX] * ipt_param_int[nd][NEQ];

  E_Int size_ssortmp = 0;
  if (ipt_param_int[0][NB_RELAX] > 1) size_ssortmp =1;
  FldArrayF ssortmp(size_tot * size_ssortmp); E_Float* ipt_ssortmp = ssortmp.begin();
  FldArrayI ssor_size(mx_sszone * threadmax_sdm * size_ssortmp); ipt_ssor_size = ssor_size.begin();

  if (ipt_param_int[0][NB_RELAX] > 1)
    {
      iptssortmp[0] = ipt_ssortmp;
      for (E_Int nd = 1; nd < nidom; nd++)
	iptssortmp[nd] = iptssortmp[ nd-1 ] + ipt_param_int[ nd-1 ][NDIMDX] * ipt_param_int[ nd-1 ][NEQ];
    }

  /// Tableau pour GMRES (on test la premiere zone 
  //
  E_Float* ipt_VectG          = NULL;
  E_Float* ipt_VectY          = NULL;
  E_Float* ipt_Hessenberg     = NULL;
  E_Float* ipt_drodmd         = NULL;
  E_Float* ipt_xmutd          = NULL;
  E_Float* ipt_norm_kry       = NULL;
  E_Float* ipt_givens         = NULL;
  E_Int num_max_vect = 2;
  E_Int size_hessenb = 1;
  E_Int ndim_xmutd   = 1;

  FldArrayF gmrestmp(size_tot); E_Float* iptgmrestmp = gmrestmp.begin();

  if (ipt_param_int[0][IMPLICITSOLVER] == 1 )
    { 
      num_max_vect = ipt_param_int[0][NB_KRYLOV];
      size_hessenb = num_max_vect*(num_max_vect-1);
      if (ipt_param_int[0][IFLOW   ] == 3) ndim_xmutd = ndimt;
      if (ipt_param_int[0][NB_RELAX] != 0)
      	{
      	  ipt_gmrestmp[0] = iptgmrestmp;
      	  for (E_Int nd = 1; nd < nidom; nd++)
      	    ipt_gmrestmp[nd] = ipt_gmrestmp[nd - 1] + ipt_param_int[nd - 1][NDIMDX] * ipt_param_int[nd - 1][NEQ];
      	}
    }
  else{ndim_drodm =1; }

  FldArrayF  VectG(num_max_vect  ); ipt_VectG = VectG.begin();
  FldArrayF  VectY(num_max_vect-1); ipt_VectY = VectY.begin();
  FldArrayF  Hessenberg(size_hessenb); ipt_Hessenberg = Hessenberg.begin();
  FldArrayF  drodmd(ndim_drodm); ipt_drodmd = drodmd.begin();
  FldArrayF   xmutd(ndim_xmutd); ipt_xmutd  = xmutd.begin();

  FldArrayF norm_kry(threadmax_sdm); ipt_norm_kry = norm_kry.begin();
  FldArrayF  givens(2 * (num_max_vect - 1)); ipt_givens = givens.begin();

  /// Tableau pour stockage senseur oscillation
  PyObject* wigArray = PyDict_GetItemString(work,"wiggle"); FldArrayF* wig;
  K_NUMPY::getFromNumpyArray(wigArray, wig, true); E_Float* iptwig = wig->begin();


  /// Tableau de travail communs explicite/implicite
  PyObject* drodmArray = PyDict_GetItemString(work,"rhs"); FldArrayF* drodm;
  K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();
  iptdrodm_transfer[0]= iptdrodm;
  for (E_Int nd = 1; nd < nidom; nd++){ iptdrodm_transfer[nd]= iptdrodm_transfer[nd-1] + ipt_param_int[nd-1][NEQ]*ipt_param_int[nd-1][NDIMDX];}
  
  // Tableau de travail coe   ( dt/vol et diags LU)
  PyObject* coeArray = PyDict_GetItemString(work,"coe"); FldArrayF* coe;
  K_NUMPY::getFromNumpyArray(coeArray, coe, true); E_Float* iptcoe = coe->begin();

  // Tableau de travail verrou omp
  PyObject* lokArray = PyDict_GetItemString(work,"verrou_omp"); FldArrayI* lok;
  K_NUMPY::getFromNumpyArray(lokArray, lok, true); E_Int* ipt_lok  = lok->begin();

  // Tableau de travail timer omp
  PyObject*  timer_omp_Array= PyDict_GetItemString(work,"TIMER_OMP"); FldArrayF* tab_timer_omp;
  K_NUMPY::getFromNumpyArray(timer_omp_Array, tab_timer_omp, true); E_Float* timer_omp = tab_timer_omp->begin();

  //
  // On recupere ipt_param_int[nd][ ILES ] infos memoire  pour vectorisation du LU (inutile et pas present dans la base en scalaire)
  //
    #ifdef E_LU_VECT
        printf("compute.cpp: Vectorisation LU pas codee."); return NULL;
    #else
        E_Int* iptludic = NULL; E_Int* iptlumax =NULL;
    #endif

  //revoir dimension nisdom_lu_max: trop gros
  E_Int mx_nidom   = mx_sszone*nidom; // nombre maximale de domaine une fois la partition Lu local active

  //printf("thread =%d\n",threadmax_sdm);
  FldArrayI ijkv_sdm(         3*threadmax_sdm); E_Int* ipt_ijkv_sdm   =  ijkv_sdm.begin();
  FldArrayI topology(         3*threadmax_sdm); E_Int* ipt_topology   =  topology.begin();
  FldArrayI ind_CL(          24*threadmax_sdm); E_Int* ipt_ind_CL     =  ind_CL.begin();
  FldArrayI ind_dm_omp(      12*threadmax_sdm); E_Int* ipt_ind_dm_omp =  ind_dm_omp.begin();

  FldArrayI tab_verrou_lhs(2*mx_nidom*threadmax_sdm); E_Int* verrou_lhs =  tab_verrou_lhs.begin();

  FldArrayF cfl( nidom*3*threadmax_sdm); E_Float* ipt_cfl = cfl.begin();

// Iteration C level
  // 
  // 
  // 
  // loop sur les iteration demander au niveau C
  // 
  // 
  //
  
   E_Int nitrun_loc = 1;
   E_Int nstep =1;
   
   E_Int flag_res=0;
   K_FASTC::souszones_list_c( ipt_param_int , ipt_param_real, ipt_ind_dm, ipt_it_lu_ssdom, iptdtloc, ipt_iskip_lu, lssiter_loc, nidom, nitrun_loc, nstep, flag_res, lexit_lu, lssiter_verif);
   flag_res=1;
   K_FASTC::souszones_list_c( ipt_param_int , ipt_param_real, ipt_ind_dm, ipt_it_lu_ssdom, iptdtloc, ipt_iskip_lu, lssiter_loc, nidom, nitrun_loc, nstep, flag_res, lexit_lu, lssiter_verif);

   E_Int display=0;
   //calcul distri si implicit ou explicit local + modulo verif
   if( (lssiter_loc ==1 || (ipt_param_int[0][EXPLOC]== 1 && ipt_param_int[0][ITYPCP]==2))  && (nitrun_loc%iptdtloc[1] == 0 || nitrun_loc == 1) )
       {
         K_FASTC::distributeThreads_c( ipt_param_int , ipt_param_real, ipt_ind_dm, nidom  , iptdtloc , mx_omp_size_int , nstep, nitrun_loc, display );
       }

   E_Int layer_mode=1;
   matvec( 
            ipt_param_int, ipt_param_real   , no_vect_test,
            nidom              , nitrun_loc       , nstep             , nssiter       , it_target, first_it,
            kimpli             , lssiter_verif    , lexit_lu          , omp_mode      , layer_mode, mpi,
            nisdom_lu_max      ,  mx_nidom        , ndimt_flt         , threadmax_sdm , mx_synchro,
            nb_pulse           ,                
            temps              ,
            ipt_ijkv_sdm       , 
            ipt_ind_dm_omp     , ipt_topology     , ipt_ind_CL        , ipt_lok, verrou_lhs, vartype, timer_omp,
            iptludic           , iptlumax         , ipt_omp          ,
            ipt_ind_dm         , ipt_it_lu_ssdom  ,
            ipt_VectG          , ipt_VectY        , iptssor           , iptssortmp    , ipt_ssor_size, ipt_drodmd, 
	    ipt_Hessenberg     , iptkrylov        , iptkrylov_transfer, ipt_norm_kry  , ipt_gmrestmp, ipt_givens,
            ipt_cfl            ,
            iptx               , ipty             , iptz              ,
            iptCellN           , iptCellN_IBC     ,
            iptro              , iptro_m1         , iptro_p1          , iptro_sfd     ,
            iptmut             , ipt_xmutd        ,
            ipti               , iptj             , iptk              , iptvol        , 
            ipti0              , iptj0            , iptk0             ,     
            ipti_df            , iptj_df          , iptk_df           , iptvol_df     , 
            iptventi           , iptventj         , iptventk          ,  
            iptrdm             ,
            iptroflt           , iptroflt2        , iptwig            , iptstat_wig   ,
            iptdrodm           , iptcoe           , iptrot            , iptdelta      , iptro_res,
            iptdrodm_transfer  ,
            ipt_param_int_tc   , ipt_param_real_tc , ipt_linelets_int, ipt_linelets_real);

  delete [] iptx; delete [] ipt_param_int;

  
  RELEASESHAREDN( wigArray       , wig  );
  RELEASESHAREDN( drodmArray     , drodm);
  RELEASESHAREDN( coeArray       , coe  );
  RELEASESHAREDN( lokArray       , lok  );
  RELEASESHAREDN( timer_omp_Array, tab_timer_omp);
  RELEASESHAREDN( dtlocArray     , dtloc);
 
  if(pyParam_int_tc  != Py_None ) { RELEASESHAREDN( pyParam_int_tc, param_int_tc);  }
  if(pyParam_real_tc != Py_None ) { RELEASESHAREDN( pyParam_real_tc, param_real_tc);}

  if (pyLinlets_int  != Py_None) { RELEASESHAREDN( pyLinlets_int , linelets_int ); }
  if (pyLinlets_real != Py_None) { RELEASESHAREDN( pyLinlets_real, linelets_real); }

  RELEASEHOOK(hook)


  Py_INCREF(Py_None);
  return Py_None;
}
