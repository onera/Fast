/*    
    Copyright 2013-2017 Onera.

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
PyObject* K_FASTS::_applyBC(PyObject* self, PyObject* args)
{
  PyObject *zone; PyObject *metric; char* var;
  if (!PyArg_ParseTuple(args, "OOs", &zone, &metric, &var )) return NULL;


  /*------*/
  /* Zone */
  /*------*/
  vector<PyArrayObject*> hook;

  /* Get numerics from zone */
  PyObject* numerics      = K_PYTREE::getNodeFromName1(zone, ".Solver#ownData");
  PyObject*          t    = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
  E_Int* ipt_param_int    = K_PYTREE::getValueAI(t, hook);
                     t    = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
  E_Float* ipt_param_real = K_PYTREE::getValueAF(t, hook);

  /*-------------------------------------*/
  /* Extraction des variables a modifier */
  /*-------------------------------------*/
  PyObject* sol_center  = K_PYTREE::getNodeFromName1(zone, "FlowSolution#Centers");
                    t   = K_PYTREE::getNodeFromName1(sol_center, var);
  E_Float*     iptro    = K_PYTREE::getValueAF(t, hook);

  /*-------------------------------------*/
  /* Extraction (x,y,z): pour forcage spatia */
  /*-------------------------------------*/
  E_Float* iptx; E_Float* ipty; E_Float* iptz;
  GET_XYZ( "GridCoordinates", zone, iptx, ipty, iptz)

  /* get metric */
  E_Float* ipti;  E_Float* iptj; E_Float* iptk; E_Float* iptvol;
  GET_TI( METRIC_TI, metric, ipt_param_int, ipti, iptj, iptk, iptvol )

  E_Float* iptventi; E_Float* iptventj; E_Float* iptventk; 
  GET_VENT( METRIC_VENT, metric, ipt_param_int, iptventi, iptventj, iptventk )
  /*----------------------------------*/
  /* Extraction des CL de la zone  */
  /*----------------------------------*/
  E_Int nd     = 0;
  E_Int lrhs   = 0;
  E_Int lcorner= 0;
  E_Int npass  = 0;
  if(lrhs == 1)  npass = 0;

  //fenetre du domaine sans les ghostcells
  FldArrayI ind_dm(6)  ; E_Int* ipt_ind_dm   = ind_dm.begin();
  ipt_ind_dm[0] = 1;
  ipt_ind_dm[1] = ipt_param_int[ IJKV];
  ipt_ind_dm[2] = 1;
  ipt_ind_dm[3] = ipt_param_int[ IJKV +1];
  ipt_ind_dm[4] = 1;
  ipt_ind_dm[5] = ipt_param_int[ IJKV +2];

  E_Int  Nbre_thread_max = 1;
#ifdef _OPENMP
  Nbre_thread_max = omp_get_max_threads();
#endif
  FldArrayI err(Nbre_thread_max);  E_Int* ierr  = err.begin();
  FldArrayI thread_topology(3*Nbre_thread_max); 
  FldArrayI   ind_dm_thread(6*Nbre_thread_max);  
  FldArrayI          ind_CL(6*Nbre_thread_max);        
  FldArrayI        shift_lu(6*Nbre_thread_max);        
  FldArrayI       ind_CL119(6*Nbre_thread_max);    
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
       E_Int ndo  = 1;

       E_Int* ipt_thread_topology = thread_topology.begin() + (ithread-1)*3;
       E_Int* ipt_ind_dm_thread   = ind_dm_thread.begin()   + (ithread-1)*6;
       E_Int* ipt_ind_CL          = ind_CL.begin()          + (ithread-1)*6;
       E_Int* ipt_ind_CL119       = ind_CL119.begin()       + (ithread-1)*6;
       E_Int* ipt_shift_lu        = shift_lu.begin( )       + (ithread-1)*6;
//       E_Float* ipt_vteta         = vteta.begin( );
//       E_Float* ipt_roteta        = roteta.begin( );

        //calcul du sous domaine a traiter par le thread 
        indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, ipt_param_int[ ITYPCP ],
                           ipt_ind_dm,
                           ipt_thread_topology, ipt_ind_dm_thread);

        ierr[ithread-1] = BCzone(nd, lrhs , lcorner, ipt_param_int, ipt_param_real, npass,
//                                 ithread, Nbre_thread_actif, ipt_thread_topology, ipt_vteta,ipt_roteta,
                                 ipt_ind_dm , ipt_ind_dm_thread, 
                                 ipt_ind_CL , ipt_ind_CL119    , ipt_shift_lu,
                                 iptro,  ipti, iptj, iptk, iptx, ipty, iptz,
                                 iptventi, iptventj, iptventk); 
  } // Fin zone // omp

    // si pointerange foireux dans BCzone, on stop
    E_Int error = ierr[0];
    for (E_Int i = 1; i < Nbre_thread_max; i++) {error = error*ierr[i];}
    if (error == 0) 
    {
      RELEASEHOOK(hook)
      PyErr_SetString(PyExc_TypeError, "applyBC: point range is invalid.");
      return NULL; 
    }
  // sortie
  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
