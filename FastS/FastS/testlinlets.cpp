/*    
    Copyright 2013-2020 Onera.

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
// test linelets : interpolation depuis le maillage de fond
//=============================================================================
PyObject* K_FASTS::_interpfromzone(PyObject* self, PyObject* args)
{

  PyObject* coordD;  
  PyObject* coordR;
  PyObject* fieldD;  
  PyObject* fieldR;
  

  E_Int nbptD; E_Int nbptR;

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "llOOOO" , &nbptD, &nbptR,&coordD, &coordR, &fieldD,&fieldR)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "iiOOOO" , &nbptD, &nbptR,&coordD, &coordR, &fieldD,&fieldR)) return NULL;
#endif

 FldArrayF* coordDf;
 E_Float* ipt_coordDf; 
 K_NUMPY::getFromNumpyArray(coordD, coordDf, true); ipt_coordDf = coordDf->begin();

 FldArrayF* coordRf;
 E_Float* ipt_coordRf; 
 K_NUMPY::getFromNumpyArray(coordR, coordRf, true); ipt_coordRf = coordRf->begin();

 FldArrayF* fieldDf;
 E_Float* ipt_fieldDf; 
 K_NUMPY::getFromNumpyArray(fieldD, fieldDf, true); ipt_fieldDf = fieldDf->begin();

 FldArrayF* fieldRf;
 E_Float* ipt_fieldRf; 
 K_NUMPY::getFromNumpyArray(fieldR, fieldRf, true); ipt_fieldRf = fieldRf->begin();

for (E_Int ir = 0; ir < nbptR; ++ir)
{
	
	E_Int imin    = 0;
    E_Int imax    = nbptD;         
    E_Int isearch = nbptD/2;

    E_Float b1;
    E_Float a1;
    E_Float alphasbeta;

    if(ipt_coordDf[nbptD-1] > ipt_coordDf[0])
    {
       while ( (imax - imin) > 1  )
       {
         if (ipt_coordRf[ir] <= ipt_coordDf[isearch])
         {
          imax = isearch;             
          isearch = imin + (imax-imin)/2;
         }
         else
         {
          imin = isearch;
          isearch = imin + (imax-imin)/2;
         }
       }
    }
    else
    {
       while ( (imax - imin) > 1  )
       {
        if (ipt_coordRf[ir] <= ipt_coordDf[isearch])
        {
         imin = isearch;
         isearch = imin + (imax-imin)/2;
       
        }
        else
        {
         imax = isearch;             
         isearch = imin + (imax-imin)/2;          
        } 
       }          
    }

    // a1 = ipt_coordRf[ir]  -ipt_coordDf[imax-1];
    // b1 = ipt_coordDf[imax]-ipt_coordDf[imax-1];

    
    // if (K_FUNC::E_abs(b1)<1.e-12) b1 = 1.e-12;
    // alphasbeta = a1/b1;

    // std::cout << "imax , coordR , coordD " << imax << " " << ipt_coordRf[ir] << " " << ipt_coordDf[imax] << " " << ipt_coordDf[imax-1] << std::endl;
    // ipt_fieldRf[ir] = (ipt_fieldDf[imax] - ipt_fieldDf[imax-1])*alphasbeta + ipt_fieldDf[imax-1];

    a1 = (ipt_fieldDf[imax] - ipt_fieldDf[imax-1])/(ipt_coordDf[imax]-ipt_coordDf[imax-1]);
    b1 = (ipt_fieldDf[imax-1]*ipt_coordDf[imax]-ipt_fieldDf[imax]*ipt_coordDf[imax-1])/(ipt_coordDf[imax]-ipt_coordDf[imax-1]);

    ipt_fieldRf[ir] = a1*ipt_coordRf[ir]+b1;

}



 Py_INCREF(Py_None);
 return Py_None;
}
