/*    
    Copyright 2013-2023 Onera.

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
#include <math.h>


# include "FastS/fastS.h"
# include "FastS/param_solver.h"
# include "string.h"
#ifdef _OPENMP
# include <omp.h>
#endif
using namespace std;
using namespace K_FLD;

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::copy_velocity_ale(PyObject* self, PyObject* args)
{
    PyObject* zones; PyObject* metrics; PyObject* work;
    E_Int omp_mode; E_Int it;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOll", &zones , &metrics, &work, &omp_mode, &it)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOii", &zones , &metrics, &work, &omp_mode, &it)) return NULL;
#endif

    /* tableau pour stocker dimension sous-domaine omp */
    E_Int threadmax_sdm = 1;
#ifdef _OPENMP
    threadmax_sdm  = omp_get_max_threads();
#endif

  //PyObject* tmp = PyDict_GetItemString(work,"MX_SYNCHRO"); E_Int mx_synchro = PyLong_AsLong(tmp); 
  //          tmp = PyDict_GetItemString(work,"MX_SSZONE");  E_Int mx_sszone  = PyLong_AsLong(tmp);

  E_Int nidom        = PyList_Size(zones);
  E_Int ndimdx       = 0;

  //printf("nombre de zone a traiter= %d\n",nidom);

  E_Float** iptventi; E_Float** iptventj; E_Float** iptventk; E_Float** iptvent_vertex;
  E_Float** iptx; E_Float** ipty; E_Float** iptz;
  E_Float** ipt_param_real; 
  E_Int** ipt_param_int;

  ipt_param_real    = new  E_Float*[nidom*8];
  iptventi          = ipt_param_real + nidom;
  iptventj          = iptventi       + nidom;
  iptventk          = iptventj       + nidom;
  iptvent_vertex    = iptventk       + nidom;
  iptx              = iptvent_vertex + nidom;
  ipty              = iptx           + nidom;
  iptz              = ipty           + nidom;
  
  ipt_param_int     = new  E_Int*[nidom];

  vector<PyArrayObject*> hook;
  vector<PyArrayObject*> motion;

  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    // check zone
    PyObject* zone    = PyList_GetItem(zones   , nd); // domaine i

    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(t, hook);
                       t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(t, hook);

    // Check metrics
    PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i

    GET_VENT( METRIC_VENT, metric, ipt_param_int[nd], iptventi[nd], iptventj[nd], iptventk[nd] )

    if(ipt_param_int[nd][ LALE ]== 0 || ipt_param_int[nd][ LALE ]== 2 || ipt_param_int[nd][ LALE ]== 3)
    { GET_XYZ( "GridCoordinates", zone, iptx[nd], ipty[nd], iptz[nd])}
    else { GET_XYZ( "GridCoordinates#Init", zone, iptx[nd], ipty[nd], iptz[nd])}

    if( ipt_param_int[nd][ NDIMDX ] > ndimdx ){ ndimdx = ipt_param_int[nd][ NDIMDX ]; } 

    PyObject* motion    = K_PYTREE::getNodeFromName1(zone   , "Motion");
                     t  = K_PYTREE::getNodeFromName1(motion , "VelocityX");
    iptvent_vertex[nd]  = K_PYTREE::getValueAF(t, hook);
  }

//
//  
//  Reservation tableau travail temporaire pour calcul du champ N+1
//

  //printf("thread =%d\n",threadmax_sdm);
  //FldArrayI compteur(     threadmax_sdm); E_Int* ipt_compteur   =  compteur.begin();
  //FldArrayI ijkv_sdm(   3*threadmax_sdm); E_Int* ipt_ijkv_sdm   =  ijkv_sdm.begin();
  FldArrayI topology(   3*threadmax_sdm); E_Int* ipt_topology   =  topology.begin();
  FldArrayI ind_dm(     6*threadmax_sdm); E_Int* ipt_ind_dm     =  ind_dm.begin();
  FldArrayI ind_dm_omp(12*threadmax_sdm); E_Int* ipt_ind_dm_omp =  ind_dm_omp.begin();

#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
# include "FastC/HPC_LAYER/INFO_SOCKET.h"

      //
      //---------------------------------------------------------------------
      // -----Boucle sur num.les domaines de la configuration
      // ---------------------------------------------------------------------
        E_Int nd_current = 0;
        for (E_Int nd = 0; nd < nidom; nd++)
          {  

            E_Int ndo   = nd; 
            
            E_Int* ipt_ind_dm_loc         = ipt_ind_dm         + (ithread-1)*6;
            ipt_ind_dm_loc[0] = 1;
            ipt_ind_dm_loc[2] = 1;
            ipt_ind_dm_loc[4] = 1;
            ipt_ind_dm_loc[1] = ipt_param_int[nd][ IJKV];
            ipt_ind_dm_loc[3] = ipt_param_int[nd][ IJKV+1];
            ipt_ind_dm_loc[5] = ipt_param_int[nd][ IJKV+2];


            if( (ipt_param_int[nd][ ITYPVENT ]== 1 && ipt_param_int[nd][ ITYPZONE ] == 1) || ipt_param_int[nd][ ITYPZONE ] == 3)
            {  ipt_ind_dm_loc[5] =1; } //champ de vitesse 2D

            E_Int* ipt_topology_socket    = ipt_topology       + (ithread-1)*3; 
            //E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm       + (ithread-1)*3; 
            E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;

             // Distribution de la sous-zone sur les threads
             //E_Int type_decoup =2;
             E_Int lmin          = 10;
             if (ipt_param_int[nd][ ITYPCP ] == 2) lmin = 4;

             indice_boucle_lu_(ndo, ithread , Nbre_thread_actif, lmin,
                              ipt_ind_dm_loc, 
                              ipt_topology_socket, ipt_ind_dm_socket );

             if (ipt_param_int[nd][ LALE ] == 2 || ipt_param_int[nd][ LALE ] == 3)
             {
               copy_ventijk_(nd, ithread, ipt_param_int[nd], ipt_param_real[nd],
                            iptx[nd], ipty[nd], iptz[nd],
                            ipt_ind_dm_loc, ipt_ind_dm_socket,
                            iptventi[nd], iptventj[nd], iptventk[nd], iptvent_vertex[nd] );
             }
            nd_current++;

          }// boucle zone 
  }  // zone OMP


  delete [] ipt_param_real;
  delete [] ipt_param_int;

  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
