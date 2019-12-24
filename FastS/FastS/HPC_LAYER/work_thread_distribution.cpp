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
# include <string.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::work_thread_distribution(PyObject* self, PyObject* args)
{
  PyObject* inddm; PyObject* topo; PyObject* ind_dm_th;
  E_Int ndo; E_Int Nthreads; E_Int lmin;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOlll", &inddm , &topo, &ind_dm_th, &ndo, &Nthreads, &lmin )) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOiii", &inddm , &topo, &ind_dm_th, &ndo, &Nthreads, &lmin )) return NULL; 
#endif


  FldArrayI*  indice_zone; 
  FldArrayI*  topo_thread; 
  FldArrayI*  indice_thread; 
  K_NUMPY::getFromNumpyArray( inddm    , indice_zone  , true); E_Int* ipt_inddm    = indice_zone->begin();
  K_NUMPY::getFromNumpyArray( topo     , topo_thread  , true); E_Int* ipt_topo     = topo_thread->begin();
  K_NUMPY::getFromNumpyArray( ind_dm_th, indice_thread, true); E_Int* ipt_inddm_th = indice_thread->begin();

  //printf("indice %d %d %d %d %d %d \n",ipt_inddm[0], ipt_inddm[1], ipt_inddm[2], ipt_inddm[3], ipt_inddm[4], ipt_inddm[5] ); 

  
  //for (E_Int ithread = 1; ithread <= Nthreads; ithread++)
  for (E_Int ithread = 1; ithread <= 1; ithread++)
     { 
            indice_boucle_lu_(ndo, ithread, Nthreads, lmin,
                              ipt_inddm,
                              ipt_topo, ipt_inddm_th );

        //printf("indicetopo %d %d %d %d \n",ipt_topo[0], ipt_topo[1], ipt_topo[2], ithread); 
        //printf("indiceLU  %d %d %d %d %d \n",ipt_inddm[0], ipt_inddm[1], ipt_inddm[2], ipt_inddm[3], ithread); 
        //printf("indiceLUth%d %d %d %d %d \n",ipt_inddm_th[0], ipt_inddm_th[1], ipt_inddm_th[2], ipt_inddm_th[3], ithread); 
     }

  RELEASESHAREDN(inddm    , indice_zone);
  RELEASESHAREDN(topo     , topo_thread);
  RELEASESHAREDN(ind_dm_th, indice_thread);


  //return Py_BuildValue("[iii]", nidom_tot, lexit_lu, lssiter_verif);
  return Py_None;
}


