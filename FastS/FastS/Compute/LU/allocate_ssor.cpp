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
 
#include "fastS.h"
#include "param_solver.h"
#include <string.h>
#ifdef _OPENMP 
#include <omp.h>
#endif
using namespace std;
using namespace K_FLD;
#include "stub.h"
extern int __activation__;

PyObject* K_FASTS::allocate_ssor(PyObject* self, PyObject* args)
{
  if (__activation__ == 0) { PyErr_SetString(PyExc_NotImplementedError, STUBMSG); return NULL; }

  PyObject *zones; PyObject *metrics; PyObject* work;
  E_Int nssiter; E_Int ompmode;

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOlOl", &zones, &metrics, &nssiter, &work, &ompmode)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOiOi", &zones, &metrics, &nssiter, &work, &ompmode)) return NULL; 
#endif
  
  PyObject* ssors  = PyList_New(0);
  PyObject* ssor; PyObject* ssortmp;
  vector<PyArrayObject*> hook;
  E_Int nidom = PyList_Size(zones), nb_subzones;
  E_Int *ipt_param_int, *ipt_ind_dm, *ipt_nidom_loc;

  PyObject* tmp = PyDict_GetItemString(work,"MX_SSZONE"); 
  E_Int mx_sszone = PyLong_AsLong(tmp);
  E_Int sizessor;
  E_Int neq;
  E_Int nfic_ij, nfic_k;

#ifdef _OPENMP
    E_Int Nbre_thread_actif = __NUMTHREADS__;
#else
    E_Int Nbre_thread_actif = 1;
#endif

  for (E_Int nd = 0; nd < nidom; nd++)
    {
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
	  nb_subzones = ipt_nidom_loc[0]; //nstep = 0

	  for (E_Int nd_subzone = 0; nd_subzone < nb_subzones; nd_subzone++)
	    {
	      if (ompmode == 1)
		{
		  E_Int       Ptomp = ipt_param_int[PT_OMP];
		  E_Int  PtrIterOmp = ipt_param_int[Ptomp];   
		  E_Int  PtZoneomp  = ipt_param_int[PtrIterOmp + nd_subzone];

		  for (E_Int i = 0; i < Nbre_thread_actif; i++)
		    {
		      if (ipt_param_int[PtZoneomp + i] != - 2)
			{
			  E_Int* ipt_ind_dm_thread = ipt_param_int + PtZoneomp +  Nbre_thread_actif + 4 + (ipt_param_int[PtZoneomp + i]) * 6;

			  sizessor += (ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] + 1 + 2 * nfic_ij) *
			    (ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] + 1 + 2 * nfic_ij) *
			    (ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] + 1 + 2 * nfic_k);
			}
		    }
		}
	      else //ompmode = 0
		{
		  E_Int* ipt_ind_dm_loc  = ipt_ind_dm + 6 * nd_subzone;
		  E_Int ipt_topology_socket[3];
		  E_Int ipt_ind_dm_thread[6];
		  E_Int lmin = 10;
		  if (ipt_param_int[ITYPCP] == 2) lmin = 4;

		  for (E_Int i = 1; i < Nbre_thread_actif + 1; i++)
		    {
		      indice_boucle_lu_(nd, i, Nbre_thread_actif, lmin,
					ipt_ind_dm_loc,
					ipt_topology_socket, ipt_ind_dm_thread);

		      if (i > ipt_topology_socket[0]*ipt_topology_socket[1]*ipt_topology_socket[2])
			break;

		      sizessor += (ipt_ind_dm_thread[1] - ipt_ind_dm_thread[0] + 1 + 2 * nfic_ij) *
			(ipt_ind_dm_thread[3] - ipt_ind_dm_thread[2] + 1 + 2 * nfic_ij) *
			(ipt_ind_dm_thread[5] - ipt_ind_dm_thread[4] + 1 + 2 * nfic_k);
		    }
		}
	    }//loop subzone
          E_Int sz_ssortmp = sizessor;
          E_Int sz_ssor    = sizessor;
          if( ipt_param_int[NB_RELAX] == 1 ) sz_ssor=1;

	  ssor    = K_NUMPY::buildNumpyArray(sz_ssor   , neq, 0, 1);  PyList_Append(ssors, ssor);
	  ssortmp = K_NUMPY::buildNumpyArray(sz_ssortmp, neq, 0, 1);  PyList_Append(ssors, ssortmp);
	  
	}// test NB_RELAX
    }//loop zone
  
  if ( nidom != 0 && (ipt_param_int[NB_RELAX] > 1 || ipt_param_int[LU_MATCH]==1 ) )
    { 
      Py_DECREF(ssor);
      Py_DECREF(ssortmp);
    }
  return ssors;
}
