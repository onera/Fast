/*    
    Copyright 2013-2024 Onera.

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
// Mettre a 1 pour un CPU timer
#define TIMER 0

# include "fastS.h"
# include "param_solver.h"
# include <string.h>
//# include <omp.h>
#if TIMER == 1
# include <ctime>
E_Float timein;
E_Float timeout;
#endif
using namespace std;
using namespace K_FLD;
#ifdef _MPI
#include <mpi.h>
#endif


#include <iostream>

//=============================================================================
// Stretching pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::_stretch(PyObject* self, PyObject* args)
{
  PyObject* coord;  
  E_Int nbpt; E_Int ityp; E_Float dx1; E_Float dx2;E_Float x2;E_Float x1;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Olffffl" , &coord, &nbpt, &x1, &x2 , &dx1, &dx2, &ityp)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oiddddi" , &coord, &nbpt, &x1, &x2 , &dx1, &dx2, &ityp)) return NULL;
#endif

 FldArrayF* coordF;
 E_Float* ipt_coordF; 
 K_NUMPY::getFromNumpyArray(coord, coordF, true); ipt_coordF = coordF->begin();

 
 // E_Float dx2 =0.0;
 // E_Int ityp =1;


 // E_Float dx2 =0.1*x2;
 // E_Float dx2 =0.6*x2;
 // E_Int ityp =2;

 E_Int inewt=1;

 k6stretch_(x1,x2,ipt_coordF,nbpt,dx1,dx2,ityp,inewt);

 Py_INCREF(Py_None);
 return Py_None;
}
