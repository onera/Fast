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
# include "fastS.h"
# include "param_solver.h"
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
PyObject* K_FASTS::_applyBC_(PyObject* self, PyObject* args)
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


  E_Int nidom    = PyList_Size(zones);

  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; 
  
  E_Float** ipt_param_real  ;
  E_Float** iptx;       E_Float** ipty;     E_Float** iptz;    E_Float** iptro;
  E_Float** ipti;       E_Float** iptj;     E_Float** iptk;    E_Float** iptvol;
  E_Float** iptventi;   E_Float** iptventj; E_Float** iptventk;

  ipt_param_int     = new E_Int*[nidom*2];
  ipt_ind_dm        = ipt_param_int   + nidom;

  iptx              = new E_Float*[nidom*12];
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
  

  /*------*/
  /* Zone */
  /*------*/
  vector<PyArrayObject*> hook;
  E_Int rk; 
  E_Int exploc; 
  E_Int* autorisation_bc = new E_Int[nidom];
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

    if (rk == 3 and exploc == 2) 
      {

	cycl = ipt_param_int[nd][NSSITER]/ipt_param_int[nd][LEVEL];
	//cout << "ccoucou" << endl;
	  
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
  FldArrayI thread_topology(3*Nbre_thread_max); 
  FldArrayI   ind_dm_thread(6*Nbre_thread_max);  
  //FldArrayI        shift_lu(6*Nbre_thread_max);        
  //FldArrayI       ind_CL119(6*Nbre_thread_max);    
  FldArrayI          ind_CL(24*Nbre_thread_max);        
//  FldArrayF       vteta(4000);    
//  FldArrayF      roteta(4000);    

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
   if( omp_mode == 0){ Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

       E_Int nd_current =0;

       for (E_Int nd = 0; nd < nidom; nd++)
          
	 {

           E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + ipt_param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
           E_Int nb_subzone     = ipt_nidom_loc [nstep -1];

           E_Int* ipt_ind_CL          = ind_CL.begin() + (ithread-1)*6;
           E_Int* ipt_ind_CL119       = ind_CL.begin() + (ithread-1)*6 + Nbre_thread_max*6 ;
           E_Int* ipt_ind_CLgmres     = ind_CL.begin() + (ithread-1)*6 + Nbre_thread_max*12;
           E_Int* ipt_shift_lu        = ind_CL.begin() + (ithread-1)*6 + Nbre_thread_max*18;
           //E_Int* ipt_ind_CL119       = ind_CL119.begin()       + (ithread-1)*6;
           //E_Int* ipt_ind_CL119       = ind_CL119.begin()       + (ithread-1)*6;
           //E_Int* ipt_shift_lu        = shift_lu.begin( )       + (ithread-1)*6;

	    //   }
	    //else {autorisation_bc = 1;}
	    //cout << "autorisation_bc : " << autorisation_bc << endl;
   
          for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
          {
            E_Int ndo   = nd;

            E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nstep-1)*6*ipt_param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

            E_Int* ipt_ind_dm_thread;
            if (omp_mode == 1)
            { 
              E_Int       Ptomp = ipt_param_int[nd][PT_OMP];
              E_Int  PtrIterOmp = ipt_param_int[nd][Ptomp +nstep -1];   
              E_Int  PtZoneomp  = ipt_param_int[nd][PtrIterOmp + nd_subzone];

              Nbre_thread_actif_loc = ipt_param_int[nd][ PtZoneomp  + Nbre_thread_actif ];
              ithread_loc           = ipt_param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
              ipt_ind_dm_thread     = ipt_param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;

              if (ithread_loc == -1) {nd_current++; continue;}
            }
            else
            { 
             E_Int* ipt_thread_topology = thread_topology.begin() + (ithread-1)*3;
                    ipt_ind_dm_thread   = ind_dm_thread.begin()   + (ithread-1)*6;

             E_Int lmin = 10;
             if (ipt_param_int[nd][ITYPCP] == 2) lmin = 4;

             indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, lmin,
                               ipt_ind_dm_loc,
                               ipt_thread_topology, ipt_ind_dm_thread);

	     //cout << ipt_ind_dm_thread[0]<<" "<<ipt_ind_dm_thread[1]<<" "<<ipt_ind_dm_thread[2]<<" "<<ipt_ind_dm_thread[3]<<endl;
            }

	    
	    if (autorisation_bc[nd] == 1)

	    {

	      //cout << "coucou zone :  " << nd << endl;
  
		ierr[ithread-1] = BCzone(nd, lrhs , nstep, lcorner, ipt_param_int[nd], ipt_param_real[nd], npass,
                                     ipt_ind_dm_loc, ipt_ind_dm_thread, 
                                     ipt_ind_CL    , ipt_ind_CL119   , ipt_ind_CLgmres, ipt_shift_lu,
                                     iptro[nd]     , ipti[nd]        , iptj[nd]       , iptk[nd]    ,
                                     iptx[nd]      , ipty[nd]        , iptz[nd]     ,
                                     iptventi[nd]  , iptventj[nd]    , iptventk[nd], iptro[nd]);


		correct_coins_(nd,  ipt_param_int[nd], ipt_ind_dm_thread , iptro[nd]);

	
	    } // autorisation

          }//loop souszone
       }//loop zone
  }//fin zone omp


  nstep = nstep_stk;

  delete [] iptx; delete [] ipt_param_int;

  RELEASEHOOK(hook)
  RELEASESHAREDN( dtlocArray  , dtloc);

  Py_INCREF(Py_None);
  return Py_None;
}
