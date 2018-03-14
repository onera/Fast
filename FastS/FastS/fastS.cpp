/*    
    Copyright 2013-2018 Onera.

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
#include "fastS.h"

int __activation__;

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyfasts [] =
{
  {"metric"              , K_FASTS::metric              ,  METH_VARARGS},
  {"_motionlaw"          , K_FASTS::_motionlaw          ,  METH_VARARGS},
  {"_movegrid"           , K_FASTS::_movegrid           ,  METH_VARARGS},
  { "cart"               , K_FASTS::cartMesh            ,  METH_VARARGS},
  {"initVars"            , K_FASTS::initVars            ,  METH_VARARGS},
  {"initNuma"            , K_FASTS::initNuma            ,  METH_VARARGS},
  {"_computePT"          , K_FASTS::_computePT          ,  METH_VARARGS},
  {"computePT_trans"     , K_FASTS::computePT_trans     ,  METH_VARARGS},
  {"_computePT_mut"      , K_FASTS::_computePT_mut      ,  METH_VARARGS},
  {"computePT_my"        , K_FASTS::computePT_my        ,  METH_VARARGS},
  {"computePT_enstrophy" , K_FASTS::computePT_enstrophy ,  METH_VARARGS},
  {"computePT_variables" , K_FASTS::computePT_variables ,  METH_VARARGS},
  {"computePT_gradient"  , K_FASTS::computePT_gradient  ,  METH_VARARGS},
  {"computePT_velocity_ale"  , K_FASTS::computePT_velocity_ale  ,  METH_VARARGS},
  {"compute_effort"      , K_FASTS::compute_effort      ,  METH_VARARGS},
  {"souszones_list"      , K_FASTS::souszones_list      ,  METH_VARARGS},
  {"stockrecup"          , K_FASTS::stockrecup          ,  METH_VARARGS},
  {"display_ss_iteration", K_FASTS::display_ss_iteration,  METH_VARARGS},
  {"_applyBC"            , K_FASTS::_applyBC            ,  METH_VARARGS},
  {"itt"                 , K_FASTS::itt                 ,  METH_VARARGS},
  {"PygetRange"          , K_FASTS::PygetRange          ,  METH_VARARGS},
  {"itt"                 , K_FASTS::itt                 ,  METH_VARARGS},
  {"compute_dpJ_dpW"       , K_FASTS::compute_dpJ_dpW       ,  METH_VARARGS},
 // {"compute_RhsIterAdjoint", K_FASTS::compute_RhsIterAdjoint,  METH_VARARGS},
 // {"rhsIter_LorURelax_Adj" , K_FASTS::LorURelaxationAdjoint ,  METH_VARARGS},
 // {"LorURelaxationAdjoint" , K_FASTS::LorURelaxationAdjoint ,  METH_VARARGS},
  {NULL, NULL}
};

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
  void initfasts();
  void initfasts()
  {
    __activation__ = K_KCORE::activation();
    Py_InitModule("fasts", Pyfasts);
    import_array();
  }
}
