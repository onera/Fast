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
 
#include "fast.h"
using namespace std;

//==============================================================================
// Check if coordinates are present
// IN: varString: varString de f
// IN: f: champ
// OUT: pos: position des coord dans varString (doit etre deja alloue a 3)
// OUT: x,y,z: pointeur sur les coords dans f
//==============================================================================
E_Int K_FAST::checkCoordinates(
  char* varString, FldArrayF& f,
  E_Int* pos,
  E_Float*& x, E_Float*& y, E_Float*& z)
{
  vector<E_Float*> fv;
  E_Int nfld = f.getNfld();
  for (E_Int i = 0; i < nfld; i++) fv.push_back(f.begin(i+1));
  return checkCoordinates(varString, fv, pos, x, y, z);
}
//==============================================================================
// Check if coordinates are present
// IN: varString: varString de f
// IN: f: champ sous forme d'un vector de ptrs
// OUT: pos: position des coord dans varString (doit etre deja alloue a 3)
// OUT: x,y,z: pointeur sur les coords dans f
//==============================================================================
E_Int K_FAST::checkCoordinates(
  char* varString, vector<E_Float*>& f,
  E_Int* pos,
  E_Float*& x, E_Float*& y, E_Float*& z)
{
  x = NULL; y = NULL; z = NULL;

  E_Int ret;
  // Verifie x
  ret = K_ARRAY::isCoordinateXPresent(varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: X coordinate is missing.");
    return 0;
  }
  x = f[ret]; pos[0] = ret+1;
  // Verifie y
  ret = K_ARRAY::isCoordinateYPresent(varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: Y coordinate is missing.");
    return 0;
  }
  y = f[ret]; pos[1] = ret+1;
  // Verifie z
  ret = K_ARRAY::isCoordinateZPresent(varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: Z coordinate is missing.");
    return 0;
  }
  z = f[ret]; pos[2] = ret+1;
  return 1;
}

//==============================================================================
// Verifie que la variable varName est presente dans f
// IN: varString: varString de f
// IN: f: champ
// OUT: pos: position de la variable dans varString
//==============================================================================
E_Int K_FAST::checkVariable(
  char* varString, FldArrayF& f,
  const char* varName, E_Int* pos, E_Float*& var)
{
  vector<E_Float*> fv;
  E_Int nfld = f.getNfld();
  for (E_Int i = 0; i < nfld; i++) fv.push_back(f.begin(i+1));
  return checkVariable(varString, fv, varName, pos, var);
}
//==============================================================================
// Verifie que la variable varName est presente dans f
// IN: varString: varString de f
// IN: f: champ sous forme de vector de ptrs
// OUT: pos: position de la variable dans varString
//==============================================================================
E_Int K_FAST::checkVariable(
  char* varString, vector<E_Float*>& f,
  const char* varName, E_Int* pos, E_Float*& var)
{
  var = NULL;
  E_Int ret;
  // Verifie var
  ret = K_ARRAY::isNamePresent(varName, varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: variable is missing.");
    return 0;
  }
  var = f[ret]; pos[0] = ret+1;
  return 1;
}

//==============================================================================
E_Int K_FAST::checkMetric(
  char* varString, FldArrayF& f,
  E_Int* pos,
  E_Float*& vol, E_Float*& surfx, E_Float*& surfy, E_Float*& surfz)
{
  vector<E_Float*> fv;
  E_Int nfld = f.getNfld();
  for (E_Int i = 0; i < nfld; i++) fv.push_back(f.begin(i+1));
  return checkMetric(varString, fv, pos, vol, surfx, surfy, surfz);
}
//==============================================================================
E_Int K_FAST::checkMetric(
  char* varString, vector<E_Float*>& f,
  E_Int* pos,
  E_Float*& vol, E_Float*& surfx, E_Float*& surfy, E_Float*& surfz)
{
  vol = NULL; surfx = NULL; surfy = NULL; surfz = NULL;

  E_Int ret;
  // Verifie vol
  ret = K_ARRAY::isNamePresent("vol", varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: vol is missing.");
    return 0;
  }
  vol = f[ret]; pos[0] = ret+1;
  // Verifie surfx
  ret = K_ARRAY::isNamePresent("surfx", varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: surfx is missing.");
    return 0;
  }
  surfx = f[ret]; pos[1] = ret+1;
  // Verifie surfy
  ret = K_ARRAY::isNamePresent("surfy", varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: surfy is missing.");
    return 0;
  }
  surfy = f[ret]; pos[2] = ret+1;
  // Verifie surfz
  ret = K_ARRAY::isNamePresent("surfz", varString);
  if (ret == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: surfz is missing.");
    return 0;
  }
  surfz = f[ret]; pos[2] = ret+1;
  return 1;
}

//=============================================================================
/*
  Verifie si les variables conservatives de l'array sont suffisantes
  pour le solveur.

  IN: varString: array var string
  IN: f: vector de la solution
  OUT: pos: position des variables dans f (commence a 1)
  OUT: ptrs sur ro, rou, rov, row, roE, cellN, sx, sy, sz (mesh velocity) 
  Retourne 0 (FAIL), 1 (SUCCESS).
*/
//=============================================================================
E_Int K_FAST::checkConsVariables(
  char* varString, FldArrayF& f,
  E_Int* pos,
  E_Float*& ro, E_Float*& rou, E_Float*& rov, E_Float*& row, E_Float*& roE,
  E_Float*& cellN, E_Float*& sx, E_Float*& sy, E_Float*& sz)
{
  vector<E_Float*> fv;
  E_Int nfld = f.getNfld();
  for (E_Int i = 0; i < nfld; i++) fv.push_back(f.begin(i+1));
  return checkConsVariables(varString, fv, pos, ro, rou, rov,
                            row, roE, cellN, sx, sy, sz);
}

//=============================================================================
E_Int K_FAST::checkConsVariables(
  char* varString, vector<E_Float*>& f,
  E_Int* pos,
  E_Float*& ro, E_Float*& rou, E_Float*& rov, E_Float*& row, E_Float*& roE,
  E_Float*& cellN, E_Float*& sx, E_Float*& sy, E_Float*& sz)
{
  ro = NULL; rou = NULL; rov = NULL; row = NULL; roE = NULL;
  cellN = NULL; sx = NULL; sy = NULL; sz = NULL;

  E_Int ret;
  // Verifie ro
  ret = K_ARRAY::isDensityPresent(varString);
  if (ret == -1) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: Density is missing.");
    return 0;
  }
  ro = f[ret]; pos[0] = ret+1;
  // Verifie rou
  ret = K_ARRAY::isMomentumXPresent(varString);
  if (ret == -1) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: MomentumX is missing.");
    return 0;
  }
  rou = f[ret]; pos[1] = ret+1;
  // Verifie rov
  ret = K_ARRAY::isMomentumYPresent(varString);
  if (ret == -1) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: MomentumY is missing.");
    return 0;
  }
  rov = f[ret]; pos[2] = ret+1;
  // Verifie row
  ret = K_ARRAY::isMomentumZPresent(varString);
  if (ret == -1) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: MomentumZ is missing.");
    return 0;
  }
  row = f[ret]; pos[3] = ret+1;
  // Verifie roE
  ret = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (ret == -1) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "compute: EnergyStagnationDensity is missing.");
    return 0;
  }
  roE = f[ret]; pos[4] = ret+1;
  // Verifie cellN
  ret = K_ARRAY::isCellNatureField1Present(varString);
  if (ret != -1) cellN = f[ret]; 

  return 1;
}
