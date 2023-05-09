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
#define K_ARRAY_UNIQUE_SYMBOL
#include "FastS/fastS.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyfasts [] =
{
  {"allocate_metric"     , K_FASTS::allocate_metric     ,  METH_VARARGS},
  {"allocate_ssor"       , K_FASTS::allocate_ssor       ,  METH_VARARGS},
  {"_movegrid"           , K_FASTS::_movegrid           ,  METH_VARARGS},
  { "cart"               , K_FASTS::cartMesh            ,  METH_VARARGS},
  {"initVars"            , K_FASTS::initVars            ,  METH_VARARGS},
  {"_matvecPT"           , K_FASTS::_matvecPT           ,  METH_VARARGS},
  {"_computePT"          , K_FASTS::_computePT          ,  METH_VARARGS},
  {"_stretch"            , K_FASTS::_stretch            ,  METH_VARARGS},
  {"_interpfromzone"     , K_FASTS::_interpfromzone     ,  METH_VARARGS},
  {"_computePT_mut"      , K_FASTS::_computePT_mut      ,  METH_VARARGS},
  {"computePT_my"        , K_FASTS::computePT_my        ,  METH_VARARGS},
  {"computePT_enstrophy" , K_FASTS::computePT_enstrophy ,  METH_VARARGS},
  {"computePT_variables" , K_FASTS::computePT_variables ,  METH_VARARGS},
  {"computePT_gradient"  , K_FASTS::computePT_gradient  ,  METH_VARARGS},
  {"computePT_velocity_ale" , K_FASTS::computePT_velocity_ale  ,  METH_VARARGS},
  {"copy_velocity_ale"   , K_FASTS::copy_velocity_ale   ,  METH_VARARGS},
  {"compute_effort"      , K_FASTS::compute_effort      ,  METH_VARARGS},
  {"stockrecup"          , K_FASTS::stockrecup          ,  METH_VARARGS},
  {"recup"               , K_FASTS::recup               ,  METH_VARARGS},
  {"recup2"              , K_FASTS::recup2              ,  METH_VARARGS},
  {"recup3"              , K_FASTS::recup3              ,  METH_VARARGS},
  {"recup3para"          , K_FASTS::recup3para          ,  METH_VARARGS},
  {"recup3para_"         , K_FASTS::recup3para_         ,  METH_VARARGS},
  {"recup3para_mpi"      , K_FASTS::recup3para_mpi      ,  METH_VARARGS},
  {"dtlocal"             , K_FASTS::dtlocal             ,  METH_VARARGS},
  {"dtlocal2"            , K_FASTS::dtlocal2            ,  METH_VARARGS},
  {"dtlocal2para"        , K_FASTS::dtlocal2para        ,  METH_VARARGS},
  {"dtlocal2para_"       , K_FASTS::dtlocal2para_       ,  METH_VARARGS},
  {"dtlocal2para_mpi"    , K_FASTS::dtlocal2para_mpi    ,  METH_VARARGS},
  {"prep_cfl"            , K_FASTS::prep_cfl            ,  METH_VARARGS},
  {"decoupe_maillage"    , K_FASTS::decoupe_maillage    ,  METH_VARARGS},
  {"display_ss_iteration", K_FASTS::display_ss_iteration,  METH_VARARGS},
  {"_applyBC"            , K_FASTS::_applyBC            ,  METH_VARARGS},
  {"itt"                 , K_FASTS::itt                 ,  METH_VARARGS},
  {"compute_dpJ_dpW"     , K_FASTS::compute_dpJ_dpW     ,  METH_VARARGS},
 
 // {"compute_RhsIterAdjoint", K_FASTS::compute_RhsIterAdjoint,  METH_VARARGS},
 // {"rhsIter_LorURelax_Adj" , K_FASTS::LorURelaxationAdjoint ,  METH_VARARGS},
 // {"LorURelaxationAdjoint" , K_FASTS::LorURelaxationAdjoint ,  METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "fasts",
        NULL,
        -1,
        Pyfasts
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_fasts();
  PyMODINIT_FUNC PyInit_fasts()
#else
  PyMODINIT_FUNC initfasts();
  PyMODINIT_FUNC initfasts()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("fasts", Pyfasts);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}
