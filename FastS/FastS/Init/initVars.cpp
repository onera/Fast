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
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "fastS.h"
//# include "converter.h"
# include "kcore.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Init fields given in nameArray to the constant value val */
//=============================================================================
PyObject* K_FASTS::initVars(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float val;
  char* varName;

#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "Osd", &array, &varName, &val))
#else
  if (!PyArg_ParseTuple(args, "Osf", &array, &varName, &val))
#endif
  {
    return NULL;
  }

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    ni, nj, nk, cn, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initVars: invalid array definition.");
    return NULL;
  }

  E_Int posvar = K_ARRAY::isNamePresent(varName, varString)+1;
  if (posvar == 0)
  {
    printf("Warning: initVars: variable name %s is not in array. Skipped...\n",varName);
  }
  else 
  { 

   E_Float* fnp = f->begin(posvar);

  E_Int ndo    =1;
  E_Int lmin   =4;
  E_Int max_thread = 1;
#ifdef _OPENMP
  max_thread = omp_get_max_threads(); // !nombre de thread maximal dans le calcul
#endif
  FldArrayI thread_topology(3*max_thread);
  FldArrayI ind_dm_omp_thread(6*max_thread);
  FldArrayI ind_dm_thread(6*max_thread);  
  E_Int  ni_loc = ni;
  E_Int  nj_loc = nj;
  //E_Int  nk_loc = nk;

  //printf("ni =%d %d %d \n",ni,nj,nk);

#pragma omp parallel default(shared)
  {
      //#ifdef E_OMP_SOUS_DOMAIN
#ifdef _OPENMP 
       E_Int  ithread           = omp_get_thread_num() +1;
       E_Int  Nbre_thread_actif = omp_get_num_threads(); // !nombre de thread actif dans cette zone
#else
       E_Int  ithread           = 1;
       E_Int  Nbre_thread_actif = 1; // !nombre de thread actif dans cette zone
#endif
      //#else
      // E_Int ithread = 1;
      // E_Int Nbre_thread_actif = 1;
      //#endif
      
      E_Int ific =2;
      E_Int kfic =2;
      if( nk == 1) kfic = 0;

      E_Int* ipt_thread_topology   = thread_topology.begin()   + 3*(ithread-1);
      E_Int* ipt_ind_dm_thread     = ind_dm_thread.begin()     + 6*(ithread-1);
      E_Int* ipt_ind_dm_omp_thread = ind_dm_omp_thread.begin() + 6*(ithread-1);

      ipt_ind_dm_thread[0] = 1;
      ipt_ind_dm_thread[2] = 1;
      ipt_ind_dm_thread[4] = 1;
      ipt_ind_dm_thread[1] = ni -2*ific-1;
      ipt_ind_dm_thread[3] = nj -2*ific-1;
      ipt_ind_dm_thread[5] = nk -2*kfic-1;
      if( nk == 1) ipt_ind_dm_thread[5] = 1;
      

      indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, lmin,
                        ipt_ind_dm_thread, 
                        ipt_thread_topology,  ipt_ind_dm_omp_thread);

     E_Int iloop1 = ipt_ind_dm_omp_thread[0];
     E_Int jloop1 = ipt_ind_dm_omp_thread[2];
     E_Int kloop1 = ipt_ind_dm_omp_thread[4];
     if( iloop1 == 1) iloop1 = iloop1 -ific;
     if( jloop1 == 1) jloop1 = jloop1 -ific;
     if( kloop1 == 1) kloop1 = kloop1 -kfic;

     E_Int iloop2 = ipt_ind_dm_omp_thread[1];
     E_Int jloop2 = ipt_ind_dm_omp_thread[3];
     E_Int kloop2 = ipt_ind_dm_omp_thread[5];
     if( iloop2 == ipt_ind_dm_thread[1]) iloop2 = iloop2 +ific+1;
     if( jloop2 == ipt_ind_dm_thread[3]) jloop2 = jloop2 +ific+1;
     if( kloop2 == ipt_ind_dm_thread[5]) kloop2 = kloop2 +kfic+1;

      for ( E_Int k = kloop1; k <= kloop2; k++)
        {
         for ( E_Int j = jloop1; j <= jloop2; j++)
           {
            for ( E_Int i = iloop1; i <= iloop2; i++)
             {

                E_Int l = (i+ific-1) + (j+ific-1)*ni_loc +(k+kfic-1)*ni_loc*nj_loc;

                fnp[l] = val;
                //printf("init %f %d%  d  %d\n", fnp[l], j,i,l);
             }
           }
        }
  }

 }
  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}
