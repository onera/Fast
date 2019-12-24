/*    
    Copyright 2013-2020 Onera.

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
# include "string.h"
#ifdef _OPENMP
# include <omp.h>
#endif
using namespace std;
using namespace K_FLD;

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::_movegrid(PyObject* self, PyObject* args)
{
  PyObject* zones; 
  if (!PyArg_ParseTuple(args, "O", &zones )) 
    return NULL;

  /* tableau pour stocker dimension sous-domaine omp */
  E_Int threadmax_sdm = 1;
#ifdef _OPENMP
  threadmax_sdm  = omp_get_max_threads();
#endif

  E_Int nidom        = PyList_Size(zones);

  //printf("nombre de zone a traiter= %d\n",nidom);

  E_Float** iptx; E_Float** ipty; E_Float** iptz;
  E_Float** iptx0; E_Float** ipty0; E_Float** iptz0;
  E_Float** ipt_param_real;

  E_Int**   ipt_param_int;

  ipt_param_real    = new  E_Float*[nidom*7];
  iptx              = ipt_param_real + nidom;
  ipty              = iptx           + nidom;
  iptz              = ipty           + nidom;
  iptx0             = iptz           + nidom;
  ipty0             = iptx0          + nidom;
  iptz0             = ipty0          + nidom;

  ipt_param_int     = new  E_Int*[nidom];

  vector<PyArrayObject*> hook;

  E_Int ndloc = 0;
  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    // check zone
    PyObject* zone    = PyList_GetItem(zones   , nd); // domaine i

    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[ndloc]  = K_PYTREE::getValueAI(t, hook);
                       t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[ndloc] = K_PYTREE::getValueAF(t, hook);

    if(ipt_param_int[ndloc][ LALE ]!= 0){ GET_XYZ( "GridCoordinates"     , zone,  iptx[ndloc],  ipty[ndloc],  iptz[ndloc]);
                                          GET_XYZ( "GridCoordinates#Init", zone, iptx0[ndloc], ipty0[ndloc], iptz0[ndloc]);
                                          ndloc = ndloc+1;
                                        }
  }
  nidom = ndloc;  // on squeeze les zone non ALE 
//
//  
//  Reservation tableau travail temporaire pour calcul du champ N+1
//

  //printf("thread =%d\n",threadmax_sdm);
  FldArrayI topology(  3*threadmax_sdm); E_Int* ipt_topology   =  topology.begin();
  FldArrayI ind_dm_omp(6*threadmax_sdm); E_Int* ipt_ind_dm_omp =  ind_dm_omp.begin();
  FldArrayI ind_dm(    6*threadmax_sdm); E_Int* ipt_ind_dm     =  ind_dm.begin();


#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
      //
      //---------------------------------------------------------------------
      // -----Boucle sur num.les domaines de la configuration
      // ---------------------------------------------------------------------
        for (E_Int nd = 0; nd < nidom; nd++)
          {  

             E_Int* ipt_topology_thread    = ipt_topology    + (ithread-1)*3; 
             E_Int* ipt_ind_dm_thread      = ipt_ind_dm      + (ithread-1)*6;
             E_Int* ind_mtr                = ipt_ind_dm_omp  + (ithread-1)*6;

             ipt_ind_dm_thread[0]=1;
             ipt_ind_dm_thread[2]=1;
             ipt_ind_dm_thread[4]=1;
             ipt_ind_dm_thread[1]= ipt_param_int[nd][IJKV  ];
             ipt_ind_dm_thread[3]= ipt_param_int[nd][IJKV+1];
             ipt_ind_dm_thread[5]= ipt_param_int[nd][IJKV+2];

             // Distribution de la sous-zone sur les threads
             E_Int lmin =4;
             indice_boucle_lu_(nd, ithread, Nbre_thread_actif, lmin,
                               ipt_ind_dm_thread, 
                               ipt_topology_thread, ind_mtr );

             if( ind_mtr[0]==ipt_ind_dm_thread[0] ) ind_mtr[0]= ind_mtr[0] -ipt_param_int[nd][ NIJK_XYZ+3]  ;
             if( ind_mtr[2]==ipt_ind_dm_thread[2] ) ind_mtr[2]= ind_mtr[2] -ipt_param_int[nd][ NIJK_XYZ+3]  ;
             if( ind_mtr[4]==ipt_ind_dm_thread[4] ) ind_mtr[4]= ind_mtr[4] -ipt_param_int[nd][ NIJK_XYZ+4]  ;
             if( ind_mtr[1]==ipt_ind_dm_thread[1] ) ind_mtr[1]= ind_mtr[1] +ipt_param_int[nd][ NIJK_XYZ+3]+1;
             if( ind_mtr[3]==ipt_ind_dm_thread[3] ) ind_mtr[3]= ind_mtr[3] +ipt_param_int[nd][ NIJK_XYZ+3]+1;
             if( ind_mtr[5]==ipt_ind_dm_thread[5] ) ind_mtr[5]= ind_mtr[5] +ipt_param_int[nd][ NIJK_XYZ+4]+1;

             move_domx_(nd,  ipt_param_int[nd], ipt_param_real[nd], ind_mtr,
                      iptx[nd]                , ipty[nd]                 , iptz[nd],  
                      iptx0[nd]               , ipty0[nd]               , iptz0[nd] );

          }// boucle zone 
  }  // zone OMP


  delete [] ipt_param_real;
  delete [] ipt_param_int;

  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
