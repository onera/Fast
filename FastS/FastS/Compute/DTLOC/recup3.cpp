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
// Idem: in place + from zone + tc compact au niveau base 
//=============================================================================
PyObject* K_FASTS::recup3(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *work;
  PyObject *pyParam_int, *pyParam_real;
  PyObject* stock;
  E_Int nstep, vartype;

  if (!PYPARSETUPLE(args,
                    "OOOOOOll", "OOOOOOii",
                    "OOOOOOll", "OOOOOOii",
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work,&stock,&vartype,&nstep))
  {
      return NULL;
  }

  vector<PyArrayObject*> hook;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);
 
   
  E_Int NoTransfert = 1;

  //// Recuperation du tableau param_int de l'arbre t et des veceurs rop et roptmp
  E_Int**   param_intt = new E_Int*[nidomR];
  E_Float** iptro_p1 = new E_Float*[nidomR];
  E_Float** iptro    = new E_Float*[nidomR];

  for (E_Int nd = 0; nd < nidomR; nd++)
     {   
       PyObject* zone = PyList_GetItem(zonesR, nd);
 
       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
       param_intt[nd]      = K_PYTREE::getValueAI(o, hook);

       o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t  = K_PYTREE::getNodeFromName1( o        , "Density_P1");
       iptro_p1[nd] = K_PYTREE::getValueAF(t, hook);

        o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t1  = K_PYTREE::getNodeFromName1( o        , "Density");
       iptro[nd]     = K_PYTREE::getValueAF(t1, hook);

     }

  /// Recuperation du tableau de stockage des valeurs
  FldArrayF* stk;
  K_NUMPY::getFromNumpyArray(stock, stk, true); E_Float* iptstk = stk->begin();

  /*-------------------------------------*/
  /* Extraction tableau int et real      */
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();
  

  E_Int cycle;

  E_Int ech  = ipt_param_int[ NoTransfert ];
  E_Int nrac = ipt_param_int[ ech +1 ];

  E_Int taille=200000000/nrac;

  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
    {

      E_Int shift_rac =  ech + 2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné

      cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];

     if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 

       {
	   E_Int donorPts[6];  E_Int idir = 2;
	   donorPts[0] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 7];
	   donorPts[1] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 8] ;
	   donorPts[2] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 9];
	   donorPts[3] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 10];
	   donorPts[4] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 11];
	   donorPts[5] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 12];
	   int dir = ipt_param_int[ech + 2 + nrac*14 + 14*irac + 13];

	   donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

	   E_Int taillefenetre;
	   taillefenetre = (donorPts[1] - donorPts[0] + 1)*(donorPts[3] - donorPts[2] + 1)*(donorPts[5] - donorPts[4] + 1); 
	  
	   if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1) 
	   //if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1) 
	     {

	       //// Recuperation de y6 en position 1 dans le tableau de stockage du raccord
	       E_Int ind=2;
	       copy_rk3local_(param_intt[NoD], donorPts , iptro[NoD], iptstk + irac*taille + 1*5*taillefenetre,ind,taillefenetre); 


	     }

       }
     else if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL]) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
       {

	   E_Int donorPts[6];  E_Int idir = 2;
	   donorPts[0] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 7];
	   donorPts[1] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 8] ;
	   donorPts[2] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 9];
	   donorPts[3] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 10];
	   donorPts[4] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 11];
	   donorPts[5] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 12];
	   int dir = ipt_param_int[ech + 2 + nrac*14 + 14*irac + 13];

	   donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

	   E_Int taillefenetre;
	   taillefenetre = (donorPts[1] - donorPts[0] + 1)*(donorPts[3] - donorPts[2] + 1)*(donorPts[5] - donorPts[4] + 1);
 
	   if (nstep%cycle==cycle/2 + cycle/4) 	   
	   //if (nstep%cycle==cycle-2) 
	
	     {
	       // cout <<"recup : "<< donorPts[0]<<" "<< donorPts[1]<<"  "<< donorPts[2]<<"  "<< donorPts[3]<<"  "<< donorPts[4]<<"  "<< donorPts[5]<<" "<<NoD<< endl;
	       //// Recuperation de y4 en position 10
	       //cout << "coucou" << endl;
	       //cout << donorPts[0] <<" "<<  donorPts[1] <<" "<<donorPts[2]<<" "<<donorPts[3] << endl;
	       E_Int ind=2;
	       copy_rk3local_(param_intt[NoD], donorPts , iptro[NoD], iptstk + irac*taille + 2*5*taillefenetre,ind,taillefenetre);

	     }







       }




    }



  //#pragma omp barrier





 RELEASESHAREDN( stock  , stk );
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
