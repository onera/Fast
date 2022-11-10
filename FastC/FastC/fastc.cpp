/*    
    Copyright 2013-2022 Onera.

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
#define K_ARRAY_UNIQUE_SYMBOL
#include "FastC/fastc.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyfastc [] =
{
  {"_motionlaw"          , K_FASTC::_motionlaw          ,  METH_VARARGS},
  {"PygetRange"          , K_FASTC::PygetRange          ,  METH_VARARGS},
  {"souszones_list"      , K_FASTC::souszones_list      ,  METH_VARARGS},
  {"distributeThreads"   , K_FASTC::distributeThreads   ,  METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "fastc",
        NULL,
        -1,
        Pyfastc
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_fastc();
  PyMODINIT_FUNC PyInit_fastc()
#else
  PyMODINIT_FUNC initfastc();
  PyMODINIT_FUNC initfastc()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("fastc", Pyfastc);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}
