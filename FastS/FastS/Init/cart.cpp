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

// Grille cartesienne reguliere
# include "fastS.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Create a cartesian mesh of nixnjxnk points 
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree. */
// ============================================================================
PyObject* K_FASTS::cartMesh(PyObject* self, PyObject* args)
{
  E_LONG ni, nj, nk;
  //E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "(ddd)(ddd)(lll)", 
                        &xo, &yo, &zo, &hi, &hj, &hk, &ni, &nj, &nk))   
#else
  if (!PyArg_ParseTuple(args, "(fff)(fff)(lll)", 
                        &xo, &yo, &zo, &hi, &hj, &hk, &ni, &nj, &nk))
#endif
  {
    return NULL;
  }
 
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cart: ni, nk, nk must be >= 1.");
    return NULL;
  }

  // Create cartesian mesh
  PyObject* tpl;
  tpl = K_ARRAY::buildArray(3, "x,y,z", ni, nj, nk);
  E_Int nijk = ni*nj*nk;

  E_Float* coord = K_ARRAY::getFieldPtr(tpl);
  E_Float* xt = coord;
  E_Float* yt = coord + nijk;
  E_Float* zt = coord + 2*nijk;

  E_Int ndo    =1;
  E_Int lmin       = 4;
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
  //printf("dx %f %f %f \n", hi,hj,hk);

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
      if( nk == 2) kfic = 0;

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

                xt[l] = xo + (i-1+ific) * hi;
                yt[l] = yo + (j-1+ific) * hj;
                zt[l] = zo + (k-1+kfic) * hk;
 
             }
           }
        }
  }

  // Return array
  return tpl;
}
