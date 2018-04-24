/*    
    Copyright 2013-2018 Onera.

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

# include "fastS.h"
# include "param_solver.h"
# include <string.h>
//# include <omp.h>
#if TIMER == 1
# include <ctime>
E_Float timein;
E_Float timeout;
#endif
using namespace std;
using namespace K_FLD;
#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::_computePT(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* work; 

  E_Int nitrun, nstep_deb, nstep_fin, layer_mode, omp_mode; 
  
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOlllllO" , &zones , &metrics, &nitrun, &nstep_deb,  &nstep_fin, &layer_mode, &omp_mode, &work)) return NULL;
#else 
  if (!PyArg_ParseTuple(args, "OOiiiiiO" , &zones , &metrics, &nitrun, &nstep_deb,  &nstep_fin, &layer_mode, &omp_mode, &work)) return NULL;
#endif

#if TIMER == 1
  clock_t c_start,c_end;
  c_start = clock();
  timein = omp_get_wtime();
#endif

  //* tableau pour stocker dimension sous-domaine omp *//
  E_Int threadmax_sdm  = __NUMTHREADS__;

  PyObject* tmp = PyDict_GetItemString(work,"MX_SSZONE");     E_Int mx_sszone     = PyLong_AsLong(tmp);  
            tmp = PyDict_GetItemString(work,"MX_SYNCHRO");    E_Int mx_synchro    = PyLong_AsLong(tmp);  
            tmp = PyDict_GetItemString(work,"FIRST_IT");      E_Int first_it      = PyLong_AsLong(tmp);
            tmp = PyDict_GetItemString(work,"mpi");           E_Int mpi           = PyLong_AsLong(tmp);


  PyObject* dtlocArray  = PyDict_GetItemString(work,"dtloc"); FldArrayI* dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc  = dtloc->begin();
  E_Int nssiter = iptdtloc[0];

  E_Int lcfl= 0;
 // 
 // 
 // 
 // partie code specifique au Transfer Data
 // 
 // 
 //
 PyObject* pyParam_int_ibc; PyObject* pyParam_real_ibc; PyObject* pyParam_int_tc;PyObject* pyParam_real_tc; PyObject* iskipArray;
 FldArrayI* param_int_ibc; FldArrayI* param_int_tc;  E_Int*  ipt_param_int_tc ; FldArrayI* iskip_lu;
 FldArrayF* param_real_tc; FldArrayF* param_real_ibc;
 E_Int*   ipt_param_int_ibc; E_Int* ipt_iskip_lu;
 E_Float* ipt_param_real_ibc; E_Float* ipt_param_real_tc;

 E_Int lssiter_loc; E_Int lexit_lu;  E_Int lssiter_verif;
 E_Int it_target = iptdtloc[4];

 if( layer_mode ==1)
 {
    lssiter_verif = 0; // par defaut, pas de calcul cfl , ni residu Newton
    if(nitrun % iptdtloc[1] == 0 || nitrun == 1) lcfl =1;

    pyParam_int_ibc = PyDict_GetItemString(work,"param_int_ibc");
    pyParam_real_ibc = PyDict_GetItemString(work,"param_real_ibc");

    K_NUMPY::getFromNumpyArray(pyParam_int_ibc, param_int_ibc, true); ipt_param_int_ibc = param_int_ibc->begin();
    K_NUMPY::getFromNumpyArray(pyParam_real_ibc, param_real_ibc, true); ipt_param_real_ibc = param_real_ibc->begin();

    pyParam_int_tc = PyDict_GetItemString(work,"param_int_tc");
    pyParam_real_tc= PyDict_GetItemString(work,"param_real_tc"); 

    if (pyParam_int_tc != Py_None)
    { K_NUMPY::getFromNumpyArray(pyParam_int_tc , param_int_tc , true); ipt_param_int_tc = param_int_tc -> begin(); }
    else{ ipt_param_int_tc = NULL; }

    if (pyParam_real_tc != Py_None)
    { K_NUMPY::getFromNumpyArray(pyParam_real_tc, param_real_tc, true); ipt_param_real_tc= param_real_tc-> begin(); }
    else{ ipt_param_real_tc = NULL; }

    if(lcfl ==1)
    {
      iskipArray = PyDict_GetItemString(work,"skip_lu");
      K_NUMPY::getFromNumpyArray(iskipArray, iskip_lu, true); ipt_iskip_lu = iskip_lu->begin();

      PyObject* tmp1 = PyDict_GetItemString(work,"lssiter_loc"); 
      if (PyLong_Check(tmp1) == true) lssiter_loc = PyLong_AsLong(tmp1);
      else lssiter_loc = PyInt_AsLong(tmp1);
    } 
 }
else
 {
  PyObject* tmp = PyDict_GetItemString(work,"lexit_lu");      lexit_lu      = PyLong_AsLong(tmp);
            tmp = PyDict_GetItemString(work,"lssiter_verif"); lssiter_verif = PyLong_AsLong(tmp);
 }
// Fin partie code specifique au Transfer Data

  E_Int nidom    = PyList_Size(zones);

  FldArrayI n0_flt(nidom); E_Int* ipt_n0_flt = n0_flt.begin();


  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; E_Int** ipt_it_lu_ssdom;
  
  E_Float** ipt_param_real  ;
  E_Float** iptx;       E_Float** ipty;     E_Float** iptz;    E_Float** iptro; E_Float** iptro_m1; E_Float** iptro_p1; E_Float** iptmut;
  E_Float** ipti;       E_Float** iptj;     E_Float** iptk;    E_Float** iptvol;
  E_Float** ipti0;      E_Float** iptj0;    E_Float** iptk0;
  E_Float** ipti_df;    E_Float** iptj_df;  E_Float** iptk_df ; E_Float** iptvol_df;
  E_Float** iptdist;    E_Float** iptventi; E_Float** iptventj; E_Float** iptventk;
  E_Float** iptrdm;
  E_Float** iptCellN  ; E_Float** ipt_cfl_zones;
  E_Float** iptro_sfd; E_Float** iptdelta; E_Float** iptfd; E_Float** iptro_zgris; E_Float** iptro_res;

  ipt_param_int     = new E_Int*[nidom*3];
  ipt_ind_dm        = ipt_param_int   + nidom;
  ipt_it_lu_ssdom   = ipt_ind_dm      + nidom;

  iptx              = new E_Float*[nidom*28];
  ipty              = iptx            + nidom;
  iptz              = ipty            + nidom;
  iptro             = iptz            + nidom;
  iptro_m1          = iptro           + nidom;
  iptro_p1          = iptro_m1        + nidom;
  iptro_sfd         = iptro_p1        + nidom;
  iptdelta          = iptro_sfd       + nidom;
  iptro_res         = iptdelta        + nidom;
  iptmut            = iptro_res       + nidom;
  ipti              = iptmut          + nidom;
  iptj              = ipti            + nidom;
  iptk              = iptj            + nidom;
  iptvol            = iptk            + nidom;
  ipti0             = iptvol          + nidom;
  iptj0             = ipti0           + nidom;
  iptk0             = iptj0           + nidom;
  ipti_df           = iptk0           + nidom;
  iptj_df           = ipti_df         + nidom;
  iptk_df           = iptj_df         + nidom;
  iptvol_df         = iptk_df         + nidom;
  iptventi          = iptvol_df       + nidom;
  iptventj          = iptventi        + nidom;
  iptventk          = iptventj        + nidom;
  iptrdm            = iptventk        + nidom;
  iptCellN          = iptrdm          + nidom; 
  ipt_param_real    = iptCellN        + nidom;
  ipt_cfl_zones     = ipt_param_real  + nidom;   //3composants: attention a la prochaine addition

  vector<PyArrayObject*> hook;

  E_Int nb_pulse      = 0;                  // =1 => pulse acoustic centree a l'origine
  E_Int ndimdx_max    = 0;
  E_Int ndimdx_sa     = 0;
  E_Int nisdom_lu_max = 0;
  E_Int ndimt         = 0;
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

    if(lcfl  == 1 && nstep_deb ==1 ) 
    { 
       PyObject* t2 = K_PYTREE::getNodeFromName1(numerics, "CFL_minmaxmoy");
       if (t2 != NULL) { ipt_cfl_zones[nd] = K_PYTREE::getValueAF(t2, hook);}
    }

    if( ipt_param_int[ nd ][ MXSSDOM_LU ] > nisdom_lu_max) nisdom_lu_max = ipt_param_int[ nd ][ MXSSDOM_LU ];

    ndimt = ndimt + ipt_param_int[nd][ NDIMDX ];

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


  /// Tableau pour stockage senseur oscillation
  PyObject* wigArray = PyDict_GetItemString(work,"wiggle"); FldArrayF* wig;
  K_NUMPY::getFromNumpyArray(wigArray, wig, true); E_Float* iptwig = wig->begin();

  /// Tableau de travail communs explicite/implicite
  PyObject* drodmArray = PyDict_GetItemString(work,"rhs"); FldArrayF* drodm;
  K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();
  
  // Tableau de travail coe   ( dt/vol et diags LU)
  PyObject* coeArray = PyDict_GetItemString(work,"coe"); FldArrayF* coe;
  K_NUMPY::getFromNumpyArray(coeArray, coe, true); E_Float* iptcoe = coe->begin();

  // Tableau de travail verrou omp
  PyObject* lokArray = PyDict_GetItemString(work,"verrou_omp"); FldArrayI* lok;
  K_NUMPY::getFromNumpyArray(lokArray, lok, true); E_Int* ipt_lok  = lok->begin();

  //
  // On recupere ipt_param_int[nd][ ILES ] infos memoire  pour vectorisation du LU (inutile et pas present dans la base en scalaire)
  //
    #ifdef E_LU_VECT
        printf("compute.cpp: Vectorisation LU pas codee."); return NULL;
    #else
        E_Int* iptludic = NULL; E_Int* iptlumax =NULL;
    #endif

  //printf("thread =%d\n",threadmax_sdm);
  FldArrayI ijkv_sdm(         3*threadmax_sdm); E_Int* ipt_ijkv_sdm   =  ijkv_sdm.begin();
  FldArrayI topology(         3*threadmax_sdm); E_Int* ipt_topology   =  topology.begin();
  FldArrayI ind_CL(           6*threadmax_sdm); E_Int* ipt_ind_CL     =  ind_CL.begin();
  FldArrayI ind_CL119(        6*threadmax_sdm); E_Int* ipt_ind_CL119  =  ind_CL119.begin();
  FldArrayI ind_dm_omp(      12*threadmax_sdm); E_Int* ipt_ind_dm_omp =  ind_dm_omp.begin();

  FldArrayF cfl( nidom*3*threadmax_sdm); E_Float* ipt_cfl    = cfl.begin();

  if(lcfl ==1 && nstep_deb == 1)
  { 
    for (E_Int i = 0;  i <  nidom*threadmax_sdm; i++){ ipt_cfl[ i*3 ] = 0;  ipt_cfl[ i*3 +1] = 100000000; ipt_cfl[ i*3 +2] = 0; }
  }

  //revoir dimension nisdom_lu_max: trop gros
  E_Int mx_nidom   = mx_sszone*nidom; // nombre maximale de domaine une fois la partition Lu local active

// SOUS PAS
  // 
  // 
  // 
  // loop sur les sous pas
  // 
  // 
  // 
  E_Int nidom_tot;
  for (E_Int nstep = nstep_deb; nstep < nstep_fin+1; ++nstep)
  {
    E_Int skip_navier = 1; //layer_mode 0: on calcul tout le temps gsdr3
    if (layer_mode ==1)
    {
       souszones_list_c( ipt_param_int , ipt_ind_dm, ipt_it_lu_ssdom, work, iptdtloc, ipt_iskip_lu, lssiter_loc, nidom, nitrun, nstep, nidom_tot, lexit_lu, lssiter_verif);

       E_Int skip = 0;
       if ( lssiter_verif == 0 && nstep == nstep_fin && ipt_param_int[0][ ITYPCP ] ==1){skip = 1;}
       skip_navier = 0;
       if (nidom_tot > 0 && skip ==0) skip_navier = 1;
    }

    //calcul Navier Stokes + appli CL
    if (skip_navier ==1)
    {

      gsdr3(ipt_param_int, ipt_param_real   ,
            nidom              , nitrun           , nstep             , nssiter       , it_target, first_it,
            kimpli             , lssiter_verif    , lexit_lu          , omp_mode      , layer_mode, mpi,
            nisdom_lu_max      ,  mx_nidom        , ndimt_flt         , threadmax_sdm , mx_synchro,
            nb_pulse           ,                
            temps              ,
            ipt_ijkv_sdm       , 
            ipt_ind_dm_omp     , ipt_topology     , ipt_ind_CL        , ipt_ind_CL119, ipt_lok,
            iptludic           , iptlumax         ,
            ipt_ind_dm         , ipt_it_lu_ssdom  ,
            ipt_cfl            ,
            iptx               , ipty             , iptz              ,
            iptCellN           ,
            iptro              , iptro_m1         , iptro_p1          , iptro_sfd     ,
            iptmut             , 
            ipti               , iptj             , iptk              , iptvol        , 
            ipti0              , iptj0            , iptk0             ,     
            ipti_df            , iptj_df          , iptk_df           , iptvol_df     , 
            iptventi           , iptventj         , iptventk          ,  
            iptrdm             ,
            iptroflt           , iptroflt2        , iptwig            , iptstat_wig   ,
            iptdrodm           , iptcoe           , iptrot            , iptdelta      , iptro_res,
            ipt_param_int_ibc       , ipt_param_real_ibc     , ipt_param_int_tc  , ipt_param_real_tc);

      if (lcfl == 1 && nstep == 1)  //mise a jour eventuelle du CFL au 1er sous-pas
      {

         for (E_Int nd = 0; nd < nidom ; nd++) 
         {
           ipt_cfl_zones[nd][0] =  ipt_cfl[0 +nd*3*threadmax_sdm];
           ipt_cfl_zones[nd][1] =  ipt_cfl[1 +nd*3*threadmax_sdm];
           ipt_cfl_zones[nd][2] =  ipt_cfl[2 +nd*3*threadmax_sdm];
           E_Int size_dom   = ipt_param_int[nd][ IJKV ]*ipt_param_int[nd][ IJKV +1]*ipt_param_int[nd][ IJKV +2];
           E_Float scale    = size_dom;
           for (E_Int ithread = 1; ithread < threadmax_sdm ; ithread++) 
              {
               E_Float* ipt_cfl_thread  = ipt_cfl + ithread*3+ nd*3*threadmax_sdm;

               ipt_cfl_zones[nd][0] = max( ipt_cfl_zones[nd][0], ipt_cfl_thread[0]);
               ipt_cfl_zones[nd][1] = min( ipt_cfl_zones[nd][1], ipt_cfl_thread[1]);
               ipt_cfl_zones[nd][2] = ipt_cfl_zones[nd][2]+ ipt_cfl_thread[2];
              }
           ipt_cfl_zones[nd][2] = ipt_cfl_zones[nd][2]/scale;
           //printf("cflmax =%f, cflmin =%f, cflmoy =%f \n",ipt_cfl_zones[0],ipt_cfl_zones[1],ipt_cfl_zones[2]);
         }
      }//test mise a jour cfl


      if((nstep == nssiter && lssiter_verif==1) || (nstep == nssiter-1 && lssiter_verif==0))  // mise a jour  du temps courant
      {
       for (E_Int nd = 0; nd < nidom ; nd++) { ipt_param_real[nd][TEMPS]= ipt_param_real[nd][TEMPS]+ ipt_param_real[nd][DTC]; }
      }
    }//test calcul NS

  }//loop nstep

  delete [] iptx; delete [] ipt_param_int;

  
  RELEASESHAREDN( wigArray    , wig  );
  RELEASESHAREDN( drodmArray  , drodm);
  RELEASESHAREDN( coeArray    , coe  );
  RELEASESHAREDN( lokArray    , lok  );
  RELEASESHAREDN( dtlocArray  , dtloc);
 
  if(layer_mode==1)
  {
    if(pyParam_int_tc  != Py_None ) { RELEASESHAREDN( pyParam_int_tc, param_int_tc);  }
    if(pyParam_real_tc != Py_None ) { RELEASESHAREDN( pyParam_real_tc, param_real_tc);}

    if(lcfl ==1) RELEASESHAREDN( iskipArray, iskip_lu);

    RELEASESHAREDN( pyParam_int_ibc , param_int_ibc);
    RELEASESHAREDN( pyParam_real_ibc, param_real_ibc);
  }

  RELEASEHOOK(hook)


#if TIMER == 1
  c_end = clock();
  E_Float duree  = (float)(c_end-c_start) / CLOCKS_PER_SEC/threadmax_sdm;
  printf("nitrun =%d, temps =%g, cpu =%g \n",nitrun_loc,temps,duree);
#endif

  Py_INCREF(Py_None);
  return Py_None;
}
