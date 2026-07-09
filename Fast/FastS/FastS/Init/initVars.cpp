/*    
    Copyright 2013-2025 Onera.

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
# include "FastS/fastS.h"
# include "FastS/param_solver.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Init fields given in nameArray to the constant value val */
//=============================================================================
PyObject* K_FASTS::initVars(PyObject* self, PyObject* args)
{
  PyObject* Pywig; PyObject* Pyparam_int; PyObject* Pydtloc;
  E_Float val;
  E_Int shift; E_Int nd_tg;

#if defined E_DOUBLEINT
#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "OOOlld" , &Pywig  , &Pyparam_int, &Pydtloc, &shift, &nd_tg,  &val )) return NULL;
#else 
  if (!PyArg_ParseTuple(args, "OOOllf" , &Pywig  , &Pyparam_int, &Pydtloc, &shift, &nd_tg, &val )) return NULL;
#endif
#else
#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "OOOiid" , &Pywig  , &Pyparam_int, &Pydtloc, &shift, &nd_tg, &val )) return NULL;
#else 
  if (!PyArg_ParseTuple(args, "OOOiif" , &Pywig  , &Pyparam_int, &Pydtloc, &shift, &nd_tg, &val )) return NULL;
#endif
#endif

  vector<PyArrayObject*> hook;
  E_Int* ipt_param_int  = K_PYTREE::getValueAI(Pyparam_int, hook);
  E_Int* iptdtloc       = K_PYTREE::getValueAI(Pydtloc, hook);

  E_Int nssiter = iptdtloc[0];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;

  E_Int nitcfg = 1;
  E_Int nbtask = ipt_omp[nitcfg-1]; 
  E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];


  FldArrayF* wig;
  K_NUMPY::getFromNumpyArray(Pywig, wig); E_Float* iptwig  = wig->begin();

#pragma omp parallel default(shared)
  {
#ifdef _OPENMP 
       E_Int  ithread           = omp_get_thread_num() +1;
       E_Int  Nbre_thread_actif = omp_get_num_threads(); 
#else
       E_Int  ithread           = 1;
       E_Int  Nbre_thread_actif = 1;
#endif
      
      E_Int ni     =ipt_param_int[NIJK  ];
      E_Int nj     =ipt_param_int[NIJK+1];
      E_Int nk     =ipt_param_int[NIJK+2];
      E_Int ific   =ipt_param_int[NIJK+3];
      E_Int kfic   =ipt_param_int[NIJK+4];
      E_Int ndimdx = ipt_param_int[NDIMDX];

        for (E_Int ntask = 0; ntask < nbtask; ntask++)
          {
             E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
             E_Int nd = ipt_omp[ pttask ];

             if(nd==nd_tg)
              {
                E_Int* ipt_inddm_omp;

                E_Int ithread_loc     = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
                ipt_inddm_omp         = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

                if (ithread_loc == -1) {continue;}
               
                E_Int iloop1 = ipt_inddm_omp[0];
                E_Int jloop1 = ipt_inddm_omp[2];
                E_Int kloop1 = ipt_inddm_omp[4];
                if( iloop1 == 1) iloop1 = iloop1 -ific;
                if( jloop1 == 1) jloop1 = jloop1 -ific;
                if( kloop1 == 1) kloop1 = kloop1 -kfic;

                E_Int iloop2 = ipt_inddm_omp[1];
                E_Int jloop2 = ipt_inddm_omp[3];
                E_Int kloop2 = ipt_inddm_omp[5];
                if( iloop2 == ipt_param_int[IJKV  ]) iloop2 = iloop2 +ific;
                if( jloop2 == ipt_param_int[IJKV+1]) jloop2 = jloop2 +ific;
                if( kloop2 == ipt_param_int[IJKV+2]) kloop2 = kloop2 +kfic;

                //printf("init %d%  d  %d %d %d %d \n",iloop1,iloop2,jloop1,jloop2,kloop1,kloop2 );
                for ( E_Int k = kloop1; k <= kloop2; k++) {
                 for ( E_Int j = jloop1; j <= jloop2; j++) {
                  for ( E_Int i = iloop1; i <= iloop2; i++) {
                    E_Int l = (i+ific-1) + (j+ific-1)*ni +(k+kfic-1)*ni*nj;

                    iptwig[l]          = val;
                    iptwig[l+ndimdx]   = val;
                    iptwig[l+ndimdx*2] = val;
                   }
                  }
                 }
              }//nd_tg
          }//task
  }//omp


  RELEASESHAREDN( Pywig       , wig  );
  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
