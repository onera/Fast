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
# include "FastS/fastS.h"
# include "FastS/param_solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
using namespace K_FLD;

//=============================================================================
/* 
   Applique les conditions aux limites d'un type donne sur la zone
   *in place* 
*/
//=============================================================================
PyObject* K_FASTS::_applyBC(PyObject* self, PyObject* args)
{
  PyObject *zones; PyObject *metrics; PyObject *work; char* var; E_Int nstep; E_Int omp_mode;  
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOlls", &zones, &metrics, &work, &nstep, &omp_mode, &var )) return NULL;
#else 
  if (!PyArg_ParseTuple(args, "OOOiis", &zones, &metrics, &work, &nstep, &omp_mode, &var )) return NULL;
#endif

  //* tableau pour stocker dimension sous-domaine omp *//
  E_Int threadmax_sdm  = __NUMTHREADS__;

  PyObject* tmp = PyDict_GetItemString(work,"MX_SSZONE");  E_Int mx_sszone = PyLong_AsLong(tmp);

  PyObject* dtlocArray  = PyDict_GetItemString(work,"dtloc"); FldArrayI* dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc  = dtloc->begin();
  E_Int nssiter = iptdtloc[0];
  E_Int* ipt_omp = iptdtloc +9 + nssiter;


  E_Int nidom    = PyList_Size(zones);

  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; 
  
  E_Float** ipt_param_real  ;
  E_Float** iptx;       E_Float** ipty;     E_Float** iptz;    E_Float** iptro;
  E_Float** ipti;       E_Float** iptj;     E_Float** iptk;    E_Float** iptvol;
  E_Float** iptventi;   E_Float** iptventj; E_Float** iptventk; E_Float** iptmut;

  ipt_param_int     = new E_Int*[nidom*2];
  ipt_ind_dm        = ipt_param_int   + nidom;

  iptx              = new E_Float*[nidom*13];
  ipty              = iptx            + nidom;
  iptz              = ipty            + nidom;
  iptro             = iptz            + nidom;
  ipti              = iptro           + nidom;
  iptj              = ipti            + nidom;
  iptk              = iptj            + nidom;
  iptvol            = iptk            + nidom;
  iptventi          = iptvol          + nidom;
  iptventj          = iptventi        + nidom;
  iptventk          = iptventj        + nidom;
  ipt_param_real    = iptventk        + nidom;
  iptmut            = ipt_param_real  + nidom;
  

  /*------*/
  /* Zone */
  /*------*/
  vector<PyArrayObject*> hook;
  E_Int rk; 
  E_Int exploc; 
  E_Int autorisation_bc[nidom];
  E_Int cycl;
  E_Int nstep_stk = nstep;

  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    // check zone
    PyObject* zone = PyList_GetItem(zones, nd); // domaine i


    /* Get numerics from zone */
    PyObject* numerics    = K_PYTREE::getNodeFromName1(zone, ".Solver#ownData");
    PyObject*          t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(t, hook);
                       t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(t, hook);

    PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i


    /*-------------------------------------*/
    /* Extraction des variables a modifier */
    /*-------------------------------------*/
    PyObject* sol_center;
    sol_center   = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
    t            = K_PYTREE::getNodeFromName1(sol_center, var);
    iptro[nd]    = K_PYTREE::getValueAF(t, hook);

    //Pointeur visqeux: mut, dist, zgris sont en acces compact
    if(ipt_param_int[nd][ IFLOW ] > 1)
      { t          = K_PYTREE::getNodeFromName1(sol_center, "ViscosityEddy");
        iptmut[nd] = K_PYTREE::getValueAF(t, hook);
      }
    else {iptmut[nd] = iptro[nd];}

    /*-------------------------------------*/
    /* Extraction (x,y,z): pour forcage spatia */
    /*-------------------------------------*/
    GET_XYZ( "GridCoordinates", zone, iptx[nd], ipty[nd], iptz[nd]) 

    /* get metric */
    E_Float* dummy;
    GET_TI(METRIC_TI  , metric, ipt_param_int[nd], ipti [nd]  , iptj [nd]  , iptk [nd]  , iptvol[nd]   )

    GET_VENT( METRIC_VENT, metric, ipt_param_int[nd], iptventi[nd], iptventj[nd], iptventk[nd] )

    ipt_ind_dm[ nd ]      =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );

    autorisation_bc[nd]=0;
    rk = ipt_param_int[0][RK];
    exploc = ipt_param_int[0][EXPLOC];

    if (exploc == 1)   // Explicit local instationnaire : on met a jour les BC en fonction du niveau en tps de la zone    
    //if (rk == 3 and exploc == 2) 
      {
	cycl = ipt_param_int[nd][NSSITER]/ipt_param_int[nd][LEVEL];
	  
	if ( nstep_stk%cycl == cycl/2 -1 and cycl != 4)
	  {
	    nstep = 1;
	    autorisation_bc[nd] = 1;
	  }
	else if (nstep_stk%cycl == cycl/2 + cycl/4 and cycl != 4)
	  {
	    nstep = 1;
	    autorisation_bc[nd] = 1;
	  }	    
	else if (nstep_stk%cycl== cycl-1 and cycl != 4 )
	  {
	    nstep = 1;
	    autorisation_bc[nd] = 1;
	  }

	else if( nstep_stk%cycl == 1 and cycl == 4 or nstep_stk%cycl == cycl/2 and cycl == 4 or nstep_stk%cycl== cycl-1 and cycl == 4 )
	  {
	    nstep = 1;
	    autorisation_bc[nd] = 1;
	  }
      } 
    else {autorisation_bc[nd] = 1;}

	    
  } // boucle zone


  /*----------------------------------*/
  /* Extraction des CL de la zone  */
  /*----------------------------------*/
  E_Int nd     = 0;
  E_Int lrhs   = 0;
  E_Int lcorner= 0;
  E_Int npass  = 0;
  if(lrhs == 1)  npass = 0;


  E_Int  Nbre_thread_max = 1;
#ifdef _OPENMP
  Nbre_thread_max = omp_get_max_threads();
#endif
  FldArrayI err(Nbre_thread_max);  E_Int* ierr  = err.begin();
  FldArrayI ind_CL(24*Nbre_thread_max); E_Int* ipt_ind_CL  = ind_CL.begin();       

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
   E_Int* ipt_ind_dm_thread; 

   E_Int nbtask = ipt_omp[nstep-1]; 
   E_Int ptiter = ipt_omp[nssiter+ nstep-1];

   for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
       E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
       E_Int nd = ipt_omp[ pttask ];

       ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
       E_Int nd_subzone      = ipt_omp[ pttask + 1 ];
       Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
       ipt_ind_dm_thread     = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

       if (ithread_loc == -1) {continue;}

       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + ipt_param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
       E_Int  nb_subzone    = ipt_nidom_loc [nstep -1];                                           //nbre sous-zone a la sousiter courante

       E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
       E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
       E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
       E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;

       E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nstep -1)*6*ipt_param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

       if (autorisation_bc[nd] == 1)
	    {

       ierr[ithread-1] = BCzone(nd, lrhs , nstep, lcorner, ipt_param_int[nd], ipt_param_real[nd], npass,
                                     ipt_ind_dm_loc, ipt_ind_dm_thread, 
                                     ipt_ind_CL_thread  , ipt_ind_CL119   , ipt_ind_CLgmres, ipt_shift_lu,
                                     iptro[nd]          , ipti[nd]        , iptj[nd]       , iptk[nd]    ,
                                     iptx[nd]           , ipty[nd]        , iptz[nd]     ,
                                     iptventi[nd]       , iptventj[nd]    , iptventk[nd], iptro[nd], iptmut[nd]);

              correct_coins_(nd,  ipt_param_int[nd], ipt_ind_dm_thread , iptro[nd]);

	    }//autorisation
     }//loop zone

  }//fin zone omp


  nstep = nstep_stk;

  delete [] iptx; delete [] ipt_param_int;

  RELEASEHOOK(hook)
  RELEASESHAREDN( dtlocArray  , dtloc);

  Py_INCREF(Py_None);
  return Py_None;
}
