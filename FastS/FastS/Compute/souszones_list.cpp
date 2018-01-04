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

# include "fastS.h"
# include "param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::souszones_list(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* work;
  E_Int nitrun; E_Int nstep;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOll", &zones , &metrics, &work, &nitrun, &nstep)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOii", &zones , &metrics, &work, &nitrun, &nstep)) return NULL; 
#endif

  E_Int lexit_lu, lssiter_verif, nidom_tot, itypcp=0;
  
  lexit_lu      = 0; // par defaut on ne skippe pas l'appel des routine LU
  lssiter_verif = 0; // par defaut, pas de calcul cfl , ni residu Newton
  nidom_tot     = 0;  

  

  // Tableau de parametre schema temporel
  PyObject* dtlocArray  = PyList_GetItem(work,5); FldArrayI*  dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc = dtloc->begin();

  vector<PyArrayObject*> hook;

  PyObject* tmp = PyList_GetItem(work,6); 
  E_Int lssiter_loc;
  if (PyLong_Check(tmp) == true) lssiter_loc = PyLong_AsLong(tmp);
  else lssiter_loc = PyInt_AsLong(tmp);

  E_Int nidom = PyList_Size(zones);
  //calcul dimension et nombre souszone
  if(nitrun % iptdtloc[1] == 0 || nitrun == 1)   //iptdtloc[1] = kmod
  {
     PyObject* iskipArray = PyList_GetItem(work,4); FldArrayI* iskip_lu;
     K_NUMPY::getFromNumpyArray(iskipArray, iskip_lu, true); E_Int* ipt_iskip_lu = iskip_lu->begin();

     E_Int nitcfg = nstep;  

     for (E_Int nd = 0; nd < nidom; nd++)
     {    
       PyObject* zone = PyList_GetItem(zones, nd);
       
       //E_Int* nijk = K_PYTREE::getValueAI(zone, hook); //taille de la zone      

       /* Get numerics from zone */
       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
       E_Int* param_int   = K_PYTREE::getValueAI(o, hook);

       //E_Int nisdom_lu_max = 0;

       PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i

      // get metric
       E_Int* ipt_ind_dm      =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );
       E_Int* ipt_it_lu_ssdom =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_ITLU) );


       E_Int ijkv_lu[3];
       ijkv_lu[0] = K_FUNC::E_max( 1, param_int[ IJKV    ]/param_int[ SIZE_SSDOM   ]);
       ijkv_lu[1] = K_FUNC::E_max( 1, param_int[ IJKV +1 ]/param_int[ SIZE_SSDOM +1]);
       ijkv_lu[2] = K_FUNC::E_max( 1, param_int[ IJKV +2 ]/param_int[ SIZE_SSDOM +2]);

       E_Int*   ipt_it_lu_ssdom_loc =  ipt_it_lu_ssdom;
       E_Int*   ipt_it_target_ssdom =  ipt_it_lu_ssdom + param_int[ MXSSDOM_LU ];
       E_Int*   ipt_it_target_old   =  ipt_it_lu_ssdom + param_int[ MXSSDOM_LU ]*2;
       E_Int*   ipt_no_lu           =  ipt_it_lu_ssdom + param_int[ MXSSDOM_LU ]*4;


       E_Int*   ipt_nisdom_residu = ipt_ind_dm + param_int[ MXSSDOM_LU ]*6*iptdtloc[0];                  //nisdom_residu(nssiter)
       E_Int*   ipt_nidom_loc     = ipt_ind_dm + param_int[ MXSSDOM_LU ]*6*iptdtloc[0] + iptdtloc[0];    //nidom_loc(nssiter)
       E_Int*   ipt_it_bloc       = ipt_ind_dm + param_int[ MXSSDOM_LU ]*6*iptdtloc[0] + iptdtloc[0]*2;  //it_bloc(nidom)


       itypcp = param_int[ ITYPCP ];
       {  
           FldArrayI ijk_lu(param_int[ MXSSDOM_LU ]*3); E_Int* ipt_ijk_lu  = ijk_lu.begin();
 
           init_ssiter_bloc_( nd                    , nitcfg               ,  iptdtloc[0] ,
                             lssiter_loc            , param_int[ ITYPCP ]  ,
                             param_int + IJKV       , ijkv_lu              , ipt_ijk_lu        , param_int + SIZE_SSDOM ,
                             param_int[ MXSSDOM_LU ], ipt_iskip_lu         ,
                             ipt_ind_dm             , ipt_nidom_loc        , ipt_it_bloc[0]    , ipt_nisdom_residu,
                             ipt_it_lu_ssdom_loc    , ipt_it_target_ssdom  , ipt_it_target_old , ipt_no_lu, param_int);
           
           nidom_tot = nidom_tot + ipt_nidom_loc[0];
       }


     } // loop zone 

     lssiter_verif = 1;

     if (nstep == iptdtloc[0] && itypcp!=2) lexit_lu  = 1;
     iptdtloc[nstep+2] = nidom_tot;

     RELEASESHAREDN(iskipArray, iskip_lu);

   }  //fin if module_verif

  else
   { nidom_tot = iptdtloc[nstep+2]; }

  RELEASESHAREDN(dtlocArray, dtloc);
  RELEASEHOOK(hook)

  return Py_BuildValue("[iii]", nidom_tot, lexit_lu, lssiter_verif);
}


