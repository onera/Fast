/*    
    Copyright 2013-2019 Onera.

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

# ifndef _FAST_FAST_H_
# define _FAST_FAST_H_

# include "kcore.h"
# include "Fortran.h"

using namespace K_FLD;

namespace K_FAST
{ 

  PyObject* _motionlaw(              PyObject* self, PyObject* args);

  //===========
  // - Check -
  //===========
  // Check a value in Numerics Dictionary
  E_Int checkNumericsValue(PyObject* numerics, const char* value,
                           E_Int& retInt, E_Float& retFloat, char*& retChar);
  /* IN: varString: chaine de variables correspondant a f
     IN: f: champ
     Rend les ptrs sur les coords et leur position dans f (commence a 1)
  */
  E_Int checkCoordinates(
    char* varString, FldArrayF& f,
    E_Int* pos,
    E_Float*& x, E_Float*& y, E_Float*& z);
  E_Int checkCoordinates(
    char* varString, std::vector<E_Float*>& f,
    E_Int* posCoords,
    E_Float*& x, E_Float*& y, E_Float*& z);
   /* IN: varString: chaine de variables correspondant a f
     IN: f: champ
     Rend un ptr sur "vol" 
  */
   E_Int checkVariable(
    char* varString, FldArrayF& f,
    const char* varName, E_Int* pos, E_Float*& var);
  E_Int checkVariable(
    char* varString, std::vector<E_Float*>& f,
    const char* varName, E_Int* pos, E_Float*& var);

  /* IN: varString: chaine de variables correspondant a f
     IN: f: champ
     Rend les ptrs sur "vol", "surfx", "surfy", "surfz" 
  */
  E_Int checkMetric(
    char* varString, FldArrayF& f,
    E_Int* pos,
    E_Float*& vol, E_Float*& surfx, E_Float*& surfy, E_Float*& surfz);
  E_Int checkMetric(
    char* varString, std::vector<E_Float*>& f,
    E_Int* pos,
    E_Float*& vol, E_Float*& surfx, E_Float*& surfy, E_Float*& surfz);

  /* IN: varString: chaine de variables correspondant a f
     IN: f: champ
     Rend des ptrs sur certains champs et leur position dans f (commence a 1)
  */
  E_Int checkConsVariables(
    char* varString, FldArrayF& f,
    E_Int* posVars,
    E_Float*& ro, E_Float*& rou, E_Float*& rov, E_Float*& row, E_Float*& roE,
    E_Float*& cellN, E_Float*& sx, E_Float*& sy, E_Float*& sz);
  E_Int checkConsVariables(
    char* varString, std::vector<E_Float*>& f,
    E_Int* posVars,
    E_Float*& ro, E_Float*& rou, E_Float*& rov, E_Float*& row, E_Float*& roE,
    E_Float*& cellN, E_Float*& sx, E_Float*& sy, E_Float*& sz);

  //==========
  // - State -
  //==========
  E_Float gamma();
  E_Float prandtl();// Prandtl number
  E_Float betaSuth(E_Float  muSuth, E_Float CSuth, E_Float TSuth);
  E_Float Cp(E_Float Cv);
  E_Float Rgp(E_Float Cv);
  
  void testFooFast();
}
#endif
