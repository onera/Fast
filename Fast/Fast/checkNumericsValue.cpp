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
 
#include "fast.h"

//=============================================================================
/*
  Verifie si la valeur existe dans le dictionnaire numerics, 
  retoure sa valeur si possible.
  Retourne 0: FAIL, 1: FOUND int, 2: FOUND float, 3: FOUND string.
*/
//=============================================================================
  E_Int K_FAST::checkNumericsValue(
    PyObject* numerics, const char* value,
    E_Int& retInt, E_Float& retFloat, char*& retChar)
{
  int ret = PyDict_Contains(numerics, PyString_FromString(value));
  if (ret <= 0) return 0;

  PyObject* o = PyDict_GetItem(numerics, PyString_FromString(value));

  if (PyString_Check(o) == true)
  {
    ret = 3;
    retChar = PyString_AsString(o);
  }
  else if (PyInt_Check(o) == true)
  {
    ret = 1;
    retInt = PyInt_AsLong(o);
  }
  else if (PyFloat_Check(o) == true)
  {
    ret = 2;
    retFloat = PyFloat_AsDouble(o);
  }
  return ret;
}
