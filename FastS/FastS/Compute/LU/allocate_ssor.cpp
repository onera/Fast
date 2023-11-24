/*    
    Copyright 2013-2024 Onera.

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
 
#include "FastS/fastS.h"
#include "FastS/param_solver.h"
#include <string.h>
#ifdef _OPENMP 
#include <omp.h>
#endif
using namespace std;
using namespace K_FLD;

PyObject* K_FASTS::allocate_ssor(PyObject* self, PyObject* args)
{
  PyObject *zones; PyObject *metrics; PyObject* work;

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOO", &zones, &metrics, &work )) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOO", &zones, &metrics, &work )) return NULL; 
#endif
  
  PyObject* ssors  = PyList_New(0);
  PyObject* ssor; PyObject* ssortmp;
  vector<PyArrayObject*> hook;
  E_Int nidom = PyList_Size(zones), nb_subzones;
  E_Int *ipt_param_int, *ipt_ind_dm, *ipt_nidom_loc;

  PyObject* dtlocArray  = PyDict_GetItemString(work,"dtloc"); FldArrayI* dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc  = dtloc->begin();
  E_Int nssiter = iptdtloc[0];
  E_Int omp_mode = iptdtloc[8];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;

  //PyObject* tmp = PyDict_GetItemString(work, "MX_SSZONE"); 
  //E_Int mx_sszone = PyLong_AsLong(tmp);
  E_Int sizessor;
  E_Int neq;
  E_Int nfic_ij, nfic_k;

#ifdef _OPENMP
    E_Int Nbre_thread_actif = __NUMTHREADS__;
#else
    E_Int Nbre_thread_actif = 1;
#endif

   E_Int nitcfg = 1;
   E_Int nbtask = ipt_omp[nitcfg-1]; 
   E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

   // loop calcul normale
   for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
      E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
      E_Int nd = ipt_omp[ pttask ];

      //ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
      E_Int nd_subzone            = ipt_omp[ pttask + 1 ];
      E_Int Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];

      PyObject* zone   = PyList_GetItem(zones  , nd);
      PyObject* metric = PyList_GetItem(metrics, nd);

      PyObject* o = K_PYTREE::getNodeFromName1(zone, ".Solver#ownData");
      o = K_PYTREE::getNodeFromName1(o, "Parameter_int"); 
      ipt_param_int = K_PYTREE::getValueAI(o, hook);

      if (ipt_param_int[NB_RELAX] > 1 || ipt_param_int[LU_MATCH]==1 )
	{
	  sizessor = 0;
	  neq = ipt_param_int[NEQ];
	  nfic_ij = ipt_param_int[NIJK + 3];
	  nfic_k  = ipt_param_int[NIJK + 4];
	  ipt_ind_dm = K_NUMPY::getNumpyPtrI(PyList_GetItem(metric, METRIC_INDM));

	  ipt_nidom_loc = ipt_ind_dm + ipt_param_int[MXSSDOM_LU] * 6 * nssiter + nssiter;

	  for (E_Int i = 0; i < Nbre_thread_actif_loc; i++)
	    {
              E_Int* ipt_ind_dm_thread   = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + i*6;

               sizessor += (ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] + 1 + 2 * nfic_ij) *
                           (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] + 1 + 2 * nfic_ij) *
                           (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] + 1 + 2 * nfic_k);
	    }
	
          E_Int sz_ssortmp = sizessor;
          E_Int sz_ssor    = sizessor;
          if (ipt_param_int[NB_RELAX] == 1) sz_ssor=1;

          //printf("sizessor %d %d %d \n",nd, sz_ssor, Nbre_thread_actif_loc);

          ssor    = K_NUMPY::buildNumpyArray(sz_ssor   , neq, 0, 1);  PyList_Append(ssors, ssor);
          ssortmp = K_NUMPY::buildNumpyArray(sz_ssortmp, neq, 0, 1);  PyList_Append(ssors, ssortmp);
       }// test NB_RELAX
    }//loop task
  
  if ( nidom != 0 && (ipt_param_int[NB_RELAX] > 1 || ipt_param_int[LU_MATCH]==1 ) )
    { 
      Py_DECREF(ssor);
      Py_DECREF(ssortmp);
    }
  return ssors;
}
