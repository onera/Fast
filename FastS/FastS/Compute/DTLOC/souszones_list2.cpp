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

#include <iostream>
# include "FastS/fastS.h"
# include "FastS/param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::souszones_list2(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* work;PyObject* zonesD;
  PyObject *pyParam_int;
  E_Int nitrun; E_Int nstep;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOOOll", &zones , &zonesD, &pyParam_int, &metrics, &work, &nitrun, &nstep)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOOOii", &zones , &zonesD, &pyParam_int, &metrics, &work, &nitrun, &nstep)) return NULL; 
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
 
           init_ssiter_bloc2_( nd                    , nitcfg               ,  iptdtloc[0] ,
                             lssiter_loc            , param_int[ ITYPCP ]  ,
                             param_int + IJKV       , ijkv_lu              , ipt_ijk_lu        , param_int + SIZE_SSDOM ,
                             param_int[ MXSSDOM_LU ], ipt_iskip_lu         ,
                             ipt_ind_dm             , ipt_nidom_loc        , ipt_it_bloc[0]    , ipt_nisdom_residu,
                             ipt_it_lu_ssdom_loc    , ipt_it_target_ssdom  , ipt_it_target_old , ipt_no_lu, param_int);
           
           nidom_tot = nidom_tot + ipt_nidom_loc[0];
       }


     } // loop zone
 
     /// Savoir si on est en explicite local ou non
     PyObject* zone = PyList_GetItem(zones, 0);
     PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
     PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "exp_local"); 
     E_Int* exploc_      = K_PYTREE::getValueAI(o, hook);
     E_Int exploc = *exploc_;

       // cout << "exploc= " << exploc << endl;

       if (exploc != 0) /// On est en explicite local
	 {
	   E_Int cycle;
	   E_Int NoTransfert = 1;

	   vector<E_Int> stockzone(nidom,0);

	   /// On recupere le tableau param_int de l'arbre t pour chaque zone
	   E_Int**   param_intt = new E_Int*[nidom];
	   for (E_Int nd = 0; nd < nidom; nd++)
	     {   
	       PyObject* zone = PyList_GetItem(zones, nd);
 	       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
	       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
	       param_intt[nd]      = K_PYTREE::getValueAI(o, hook);
	     }

	   E_Int nssiter = param_intt[0][NSSITER];

	   /// On recupère param_int de l'arbre tc
	   FldArrayI* param_int;
	   E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
	   E_Int* ipt_param_int = param_int->begin();

	   E_Int ech  = ipt_param_int[ NoTransfert ];
	   E_Int nrac = ipt_param_int[ ech +1 ];


	   for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
	     {

	       E_Int shift_rac =  ech + 2 + irac;

	       E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	       E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné

	       cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];

	       if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL] and (nitcfg%cycle)==cycle/4)  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse 
		 {

		   stockzone[NoD]=stockzone[NoD]+1;
		  
		   PyObject* metric = PyList_GetItem(metrics, NoD); 
		   E_Int* ipt_ind_dm      =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );

		   //cout << "param_int[NoD][ MXSSDOM_LU ]= " << param_intt[NoD][ MXSSDOM_LU ] << endl;

		   E_Int* ipt_ind_dm_loc  = ipt_ind_dm  + (nitcfg-1)*6*param_intt[NoD][ MXSSDOM_LU ] + 6*(stockzone[NoD]-1);
		   ipt_ind_dm_loc[0] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 7];
		   ipt_ind_dm_loc[1] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 8] ;
		   ipt_ind_dm_loc[2] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 9];
		   ipt_ind_dm_loc[3] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 10];
		   ipt_ind_dm_loc[4] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 11];
		   ipt_ind_dm_loc[5] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 12];
		   int dir = ipt_param_int[ech + 2 + nrac*14 + 14*irac + 13];
		   ipt_ind_dm_loc[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = ipt_ind_dm_loc[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

		   E_Int* ipt_nidom_loc = ipt_ind_dm + param_intt[NoD][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
		   ipt_nidom_loc [nitcfg-1] =stockzone[NoD] ;                                          		  
		   

		 }
		   




		 }


	     }




	 







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


