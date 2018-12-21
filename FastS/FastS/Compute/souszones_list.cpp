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

# include "fastS.h"
# include "param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

// 
// 
// Souszonelist (C layer)
// 
// 
void K_FASTS::souszones_list_c( E_Int**& param_int, E_Int**& ipt_ind_dm, E_Int**& ipt_it_lu_ssdom, PyObject* work,
                                E_Int* iptdtloc   , E_Int* ipt_iskip_lu, E_Int lssiter_loc       , E_Int nidom      , 
                                E_Int nitrun      , E_Int nstep        , E_Int& nidom_tot        , E_Int& lexit_lu  , E_Int& lssiter_verif )
{

  lexit_lu      = 0; // par defaut on ne skippe pas l'appel des routine LU
  lssiter_verif = 0; // par defaut, pas de calcul cfl , ni residu Newton
  nidom_tot     = 0;  

  vector<PyArrayObject*> hook;

  //calcul dimension et nombre souszone
  if(nitrun % iptdtloc[1] == 0 || nitrun == 1)   //iptdtloc[1] = kmod
  {
     E_Int itypcp;
     for (E_Int nd = 0; nd < nidom; nd++)
     {    
       E_Int ijkv_lu[3];
       ijkv_lu[0] = K_FUNC::E_max( 1, param_int[nd][ IJKV    ]/param_int[nd][ SIZE_SSDOM   ]);
       ijkv_lu[1] = K_FUNC::E_max( 1, param_int[nd][ IJKV +1 ]/param_int[nd][ SIZE_SSDOM +1]);
       ijkv_lu[2] = K_FUNC::E_max( 1, param_int[nd][ IJKV +2 ]/param_int[nd][ SIZE_SSDOM +2]);

       E_Int*   ipt_it_lu_ssdom_loc =  ipt_it_lu_ssdom[nd];
       E_Int*   ipt_it_target_ssdom =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ];
       E_Int*   ipt_it_target_old   =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ]*2;
       E_Int*   ipt_no_lu           =  ipt_it_lu_ssdom[nd] + param_int[nd][ MXSSDOM_LU ]*4;


       E_Int*   ipt_nisdom_residu = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*iptdtloc[0];                  //nisdom_residu(nssiter)
       E_Int*   ipt_nidom_loc     = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*iptdtloc[0] + iptdtloc[0];    //nidom_loc(nssiter)
       E_Int*   ipt_it_bloc       = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*iptdtloc[0] + iptdtloc[0]*2;  //it_bloc(nidom)


       itypcp = param_int[nd][ ITYPCP ];
       {  
           FldArrayI ijk_lu(param_int[nd][ MXSSDOM_LU ]*3); E_Int* ipt_ijk_lu  = ijk_lu.begin();
 
           init_ssiter_bloc_( nd                        , nstep                    ,  iptdtloc[0] ,
                             lssiter_loc                , param_int[nd][ ITYPCP ]  ,
                             param_int[nd] + IJKV       , ijkv_lu                  , ipt_ijk_lu    , param_int[nd] + SIZE_SSDOM ,
                             param_int[nd][ MXSSDOM_LU ], ipt_iskip_lu         ,
                             ipt_ind_dm[nd]             , ipt_nidom_loc        , ipt_it_bloc[0]    , ipt_nisdom_residu,
                             ipt_it_lu_ssdom_loc        , ipt_it_target_ssdom  , ipt_it_target_old , ipt_no_lu, param_int[nd]);
           
           nidom_tot = nidom_tot + ipt_nidom_loc[0];
       }


     } // loop zone 

     lssiter_verif = 1;

     if (nstep == iptdtloc[0] && itypcp!=2) lexit_lu  = 1;
     iptdtloc[nstep+5] = nidom_tot;

   }  //fin if module_verif

  else
   { nidom_tot = iptdtloc[nstep+5]; }

  RELEASEHOOK(hook)

  return;
}

// -----------------------------------------------------------------------------------
// 
// Souszonelist (python layer)
// 
// -----------------------------------------------------------------------------------
PyObject* K_FASTS::souszones_list(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* work;
  E_Int nitrun; E_Int nstep; E_Int distrib_omp;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOlll", &zones , &metrics, &work, &nitrun, &nstep, &distrib_omp)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOiii", &zones , &metrics, &work, &nitrun, &nstep, &distrib_omp)) return NULL; 
#endif
  
  E_Int lexit_lu, nidom_tot, lssiter_verif;

  // Tableau de parametre schema temporel
  PyObject* dtlocArray  = PyDict_GetItemString(work, "dtloc"); FldArrayI*  dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc = dtloc->begin();

  E_Int nidom = PyList_Size(zones);

  PyObject* tmp = PyDict_GetItemString(work,"lssiter_loc"); 
  E_Int lssiter_loc;
  if (PyLong_Check(tmp) == true) lssiter_loc = PyLong_AsLong(tmp);
  else lssiter_loc = PyInt_AsLong(tmp);

  tmp = PyDict_GetItemString(work,"MX_OMP_SIZE_INT"); 
  E_Int mx_omp_size_int;
  if (PyLong_Check(tmp) == true) mx_omp_size_int = PyLong_AsLong(tmp);
  else mx_omp_size_int = PyInt_AsLong(tmp);

  PyObject* iskipArray = PyDict_GetItemString(work,"skip_lu"); FldArrayI* iskip_lu;
  K_NUMPY::getFromNumpyArray(iskipArray, iskip_lu, true); E_Int* ipt_iskip_lu = iskip_lu->begin();

 
  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; E_Int** ipt_it_lu_ssdom;
  ipt_param_int     = new E_Int*[nidom*3];
  ipt_ind_dm        = ipt_param_int   + nidom;
  ipt_it_lu_ssdom   = ipt_ind_dm      + nidom;

  vector<PyArrayObject*> hook;

  for (E_Int nd = 0; nd < nidom; nd++)
  {    
    PyObject* zone = PyList_GetItem(zones, nd);
       
    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);


    // get metric
    PyObject* metric     = PyList_GetItem(metrics, nd); // metric du domaine i
    ipt_ind_dm[nd]       =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );
    ipt_it_lu_ssdom[nd]  =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_ITLU) );
  }

  souszones_list_c( ipt_param_int , ipt_ind_dm, ipt_it_lu_ssdom, work, iptdtloc, ipt_iskip_lu, lssiter_loc, nidom, nitrun, nstep, nidom_tot, lexit_lu, lssiter_verif);

  //calcul distri si implicit ou explicit local + modulo verif
  if( distrib_omp==1 || (lssiter_loc ==1 || (ipt_param_int[0][EXPLOC]== 1 && ipt_param_int[0][ITYPCP]==2))  && (nitrun%iptdtloc[1] == 0 || nitrun == 1) )
  {
    E_Int display = 0;
    //if (nstep==1) display = 1;
    distributeThreads_c( ipt_param_int , ipt_ind_dm, nidom  , iptdtloc[0] , mx_omp_size_int , nstep, nitrun, display );
  }

  PyObject* dico = PyDict_New();
  tmp = Py_BuildValue("i", nidom_tot);
  PyDict_SetItemString(dico , "nidom_tot"    , tmp);

  tmp = Py_BuildValue("i", lexit_lu);
  PyDict_SetItemString(dico , "lexit_lu"     , tmp);

  tmp = Py_BuildValue("i", lssiter_verif);
  PyDict_SetItemString(dico , "lssiter_verif", tmp);

  delete [] ipt_param_int;

  RELEASESHAREDN(dtlocArray, dtloc);
  RELEASESHAREDN( iskipArray, iskip_lu);
  RELEASEHOOK(hook)

  return dico;
}
