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
# include "fastS.h"
# include "param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Fonction de recup des valeurs apres transferts
//=============================================================================
PyObject* K_FASTS::recup(PyObject* self, PyObject* args)
{

  PyObject* zones; PyObject* stock; E_Int nstep; 
  E_Int taille_tabs;
 
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOll", &zones , &stock, &nstep, &taille_tabs)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOii", &zones , &stock, &nstep, &taille_tabs)) return NULL; 
#endif

  vector<PyArrayObject*> hook;

  E_Int nidom = PyList_Size(zones);

  E_Int**   param_int   = new E_Int*[nidom];
  E_Float** iptro       = new E_Float*[nidom];
  E_Float** iptro_p1       = new E_Float*[nidom];

   
  /// Recuperation du tableau de stockage des valeurs
  FldArrayF* stk;
  K_NUMPY::getFromNumpyArray(stock, stk, true); E_Float* iptstk = stk->begin();

   
  /// Recuperation de param_int
  for (E_Int nd = 0; nd < nidom; nd++)
     {   
       PyObject* zone = PyList_GetItem(zones, nd);
 
       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
       param_int[nd]      = K_PYTREE::getValueAI(o, hook);

        o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t1  = K_PYTREE::getNodeFromName1( o        , "Density");
       iptro[nd]     = K_PYTREE::getValueAF(t1, hook);

        o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t2  = K_PYTREE::getNodeFromName1( o        , "Density_P1");
       iptro_p1[nd]  = K_PYTREE::getValueAF(t2, hook);

     }


   E_Int taille=taille_tabs/nidom;
   E_Int cycle;
   
     
   if(param_int[0][RK]==3 and param_int[0][EXPLOC]==2)
     {
       for (E_Int nd = 0; nd < nidom; nd++)
	 {
	   cycle=param_int[0][NSSITER]/param_int[nd][LEVEL];
	  
	
	   if (cycle==4) //////////////////////////////////// Zones de plus petit niveau en temps ////////////////////////////////////////////////// 
	     {
	       
	       if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1 ) /// Recuperation de y2 ou y6
		 {
		   if(param_int[nd][LEVELG]==0) /// savoir si on stocke les colonnes de gauche ou de droite
		     {
	   
		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;E_Int pos=1;		       		  
		       /// recup de y6 en position 1 dans le tableau de stockage
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  


		     }
		   else if(param_int[nd][LEVELD]==0) /// savoir quelles colonnes on prend
		     {
		       //cout << nd << endl;
		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;E_Int pos=1;		       		  
		       /// recup de y6 en position 1 dans le tableau de stockage
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);
		   		   
		     }	     


		   else
		     {
		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;E_Int pos=2;		       		  
		       /// recup de y6 en position 2 dans le tableau de stockage
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille); 

		       pos=3;
		       /// recup de y6 en position 3 dans le tableau de stockage
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);
		     }


		 }    
	     }
	  	  	 
      else if (cycle == param_int[nd][NSSITER]) /// Zones de + gd pas de temps 
	     {

	       if (param_int[nd][LEVELG]==0) /// La zone est situee sur le bord gauche du domaine et transmet des infos à droite
		 {
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];		   
		    
		   if (nstep%cycle == cycle-2)
		     {
		       /// Récupération de y4 en position 4
		       E_Int ind=2;E_Int pos=4;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);
		     }

		 }
	       if (param_int[nd][LEVELD]==0) /// La zone est situee sur le bord droit du domaine et transmet des infos à gauche
		 {
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  2;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];		   
		    
		   if (nstep%cycle == cycle-2)
		     {
		       /// Récupération de y4 en position 4
		       E_Int ind=2;E_Int pos=4;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);
		     }

		 }

	     }
	   
       else
	 {

	       if (param_int[nd][LEVELG] < param_int[nd][LEVEL] and param_int[nd][LEVEL] < param_int[nd][LEVELD])

		 {

		       E_Int idir = 2; E_Int ind_loop[6];
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];

		       if (nstep%cycle==cycle-2) 
			 {

		           //// Recuperation de y4 en position 10
			   E_Int ind=2;E_Int pos=5;
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

			   if ((nstep/cycle)%2==1)

			     {
			       ind_loop[0]=  1;
			       ind_loop[1]=  2;

			       //// Recuperation de y6 en position 12
			       ind=2; pos=6;
			       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

			     }


			 }


		 }


	       if (param_int[nd][LEVELG] > param_int[nd][LEVEL] and param_int[nd][LEVEL] > param_int[nd][LEVELD])

		 {

		       E_Int idir = 2; E_Int ind_loop[6];
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];


		       if (nstep%cycle==cycle-2) 
			 {

		           //// Recuperation de y4 en position 10
			   E_Int ind=2;E_Int pos=5;
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

			   if ((nstep/cycle)%2==1)

			     {
			       ind_loop[0]=  param_int[nd][IJKV]-1;
			       ind_loop[1]=  param_int[nd][IJKV];

			       //// Recuperation de y6 en position 12
			       ind=2; pos=6;
			       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

			     }


			 }



		 }






	 }
	  	 
	}
     }

  
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


     if (param_int[0][RK]==2 and param_int[0][EXPLOC]==4) // RK2 local test
      {      
      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  cycle=param_int[0][NSSITER]/param_int[nd][LEVEL];

	  if ((nstep+1)%cycle==0 and param_int[nd][LEVEL] < cycle) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
  	   
		  E_Int ind_loop[6]; E_Int ind=2; E_Int idir=2;		  
		   /// Recuperation des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-5;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

    
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  E_Int ind_loop[6]; E_Int ind=2; E_Int idir=2;
		   /// Recuperation des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  6;
		   ind_loop[2]=  1;//- param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

		   E_Int pos=1;       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);

		   
		}
	    }

	  
	  if ((nstep+1)%cycle==0 and param_int[nd][LEVEL] == cycle and (nstep/cycle)%2==1) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      //cout << "coucou2" << endl;

	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// Recup des valeurs
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;E_Int pos=8;
		   ind_loop[0]=  1;
		   ind_loop[1]=  2;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);

   
		  //// Recup des valeurs
		   pos=9;
		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];

		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);

   
		}
	    }


	  	  

	}
      }


 delete [] param_int;
  RELEASESHAREDN( stock  , stk );
 
  Py_INCREF(Py_None);
  return Py_None;  
} 
