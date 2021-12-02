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
 
#include "fastc.h"

//=============================================================================
/*
  Verifie si la valeur existe dans le dictionnaire numerics, 
  retoure sa valeur si possible.
  Retourne 0: FAIL, 1: FOUND int, 2: FOUND float, 3: FOUND string.
*/
//=============================================================================
  E_Int K_FASTC::checkNumericsValue(
    PyObject* numerics, const char* value,
    E_Int& retInt, E_Float& retFloat, char*& retChar)
{
  PyObject* val = PyString_FromString(value);
  int ret = PyDict_Contains(numerics, val);
  if (ret <= 0) return 0;

  PyObject* o = PyDict_GetItem(numerics, val);
  if (PyString_Check(o))
  {
    ret = 3;
    retChar = PyString_AsString(o);
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(o))
  {
    ret = 3;
    retChar = (char*)PyUnicode_AsUTF8(o);
  }
#endif
  else if (PyInt_Check(o))
  {
    ret = 1;
    retInt = PyInt_AsLong(o);
  }
  else if (PyFloat_Check(o))
  {
    ret = 2;
    retFloat = PyFloat_AsDouble(o);
  }
  return ret;
}
