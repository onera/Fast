/*    
    Copyright 2013-2016 Onera.

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
#ifdef vtune
//#include "ittnotify.h"
#endif

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Enable start and pause of data collection for Intel profiling tools. */
// ============================================================================
PyObject* K_FASTS::itt(PyObject* self, PyObject* args)
{

  E_Int flag;

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "l", &flag)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "i", &flag)) return NULL;
#endif

  //printf( "itt");
#ifdef vtune
//  if(flag==0) __itt_resume();
//  else        __itt_pause();
#endif

  Py_INCREF(Py_None);
  return Py_None;
}
