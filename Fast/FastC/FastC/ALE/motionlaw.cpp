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


# include "FastC/fastc.h"
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
PyObject* K_FASTC::_motionlaw(PyObject* self, PyObject* args)
{
  PyObject* zones; E_Float teta;  E_Float tetap;

#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "Odd", &zones, &teta, &tetap)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Off", &zones, &teta, &tetap)) return NULL;
#endif

  /* tableau pour stocker dimension sous-domaine omp */
  E_Int nidom        = PyList_Size(zones);

  E_Float** ipt_param_real    = new  E_Float*[nidom];
  E_Int**   ipt_param_int     = new  E_Int*[nidom];

  vector<PyArrayObject*> hook;

  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    // check zone
    PyObject* zone    = PyList_GetItem(zones   , nd); // domaine i

    /* Get numerics from zone */
    PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*       t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]  = K_PYTREE::getValueAI(t, hook);
                    t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd] = K_PYTREE::getValueAF(t, hook);

    if(ipt_param_int[nd][LALE]==1){ ipt_param_real[nd][ROT_TETA]=teta; ipt_param_real[nd][ROT_TETAP]=tetap;}
    
  }

  delete [] ipt_param_real;
  delete [] ipt_param_int;

  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
