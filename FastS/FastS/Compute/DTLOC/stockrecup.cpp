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
// Fonction de stockage des flux
//=============================================================================
PyObject* K_FASTS::stockrecup(PyObject* self, PyObject* args)
{

  PyObject* zones; PyObject* stock; PyObject* work;
  E_Int nstep; PyObject* drodmstock;

 
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOOl", &zones , &stock, &drodmstock, &work, &nstep)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOOi", &zones , &stock, &drodmstock, &work,  &nstep)) return NULL; 
#endif

  vector<PyArrayObject*> hook;

  E_Int nidom = PyList_Size(zones);

  E_Int**   param_int = new E_Int*[nidom];
  E_Float** param_real= new E_Float*[nidom];
  E_Float** iptro_p1 = new E_Float*[nidom];
  E_Float** iptro    = new E_Float*[nidom];

  
  /// Recuperation du tableau de stockage des valeurs
  FldArrayF* stk;
  K_NUMPY::getFromNumpyArray(stock, stk, true); E_Float* iptstk = stk->begin();

  /// Recuperation du tableau de stockage des flux
  FldArrayF* drodmstk;
  K_NUMPY::getFromNumpyArray(drodmstock,drodmstk, true); E_Float* iptdrodmstk = drodmstk->begin();

  /// Recuperation du tableau des flux (drodm)
  PyObject* drodmArray = PyList_GetItem(work,2); FldArrayF* drodm;
  K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();
   
  /// Recuperation de param_int, param_real,rop et rop_tmp
  for (E_Int nd = 0; nd < nidom; nd++)
     {   
       PyObject* zone = PyList_GetItem(zones, nd);
 
       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
       param_int[nd]      = K_PYTREE::getValueAI(o, hook);

	               o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
       param_real[nd]     = K_PYTREE::getValueAF(o, hook);

       o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t  = K_PYTREE::getNodeFromName1( o        , "Density_P1");
       iptro_p1[nd] = K_PYTREE::getValueAF(t, hook);

        o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t1  = K_PYTREE::getNodeFromName1( o        , "Density");
       iptro[nd]     = K_PYTREE::getValueAF(t1, hook);
      }  

  E_Int a;
  a=0;
  E_Int shift_zone[nidom];

     for (E_Int nd = 0; nd < nidom; nd++)
       {
	 shift_zone[nd]=a;
	 a=a+param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];	 
       }



      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  if (nstep%(param_int[0][NSSITER]/param_int[nd][LEVEL])==1) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	      
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //cout <<ipt_param_int[nd][IJKV]-2 <<endl;
		  //cout <<ipt_param_int[nd][IJKV]   <<endl;
		  //cout <<ipt_param_int[nd][IJKV+1] <<endl;
		  //cout <<"nzone = "<<nd<<"stockage"<<endl;
		  
		  //// stockage des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  param_int[nd][IJKV]-2;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop,iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);

		  
		   /// stockage des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-3;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1]+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   /// Interpolation
		   interpolation_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);
		   //cout<<"nd= "<<nd<<"  "<<"nstep= "<<nstep<<endl;
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// stockage des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  1;
		   ind_loop[1]=  3;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);

		   /// stockage des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  4;
		   ind_loop[2]=  1- param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1]+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   /// Interpolation
		   interpolation_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);
		     }	     
	    }    

	  if (((nstep+1)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //cout <<ipt_param_int[nd][IJKV]-2 <<endl;
		  //cout <<ipt_param_int[nd][IJKV]   <<endl;
		  //cout <<ipt_param_int[nd][IJKV+1] <<endl;
		  //cout <<ipt_param_int[nd][IJKV+2] <<endl;
		  //cout <<"nzone = "<<nd<<"recuperation"<<endl;
		  
		  //// Recuperation des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  param_int[nd][IJKV]-2;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);

		  
		   /// Recuperation des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-3;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1]+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);
		   //cout<<"nd= "<<nd<<"  "<<"nstep= "<<nstep<<endl;
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// Recuperation des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  1;
		   ind_loop[1]=  3;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);

		   /// Recuperation des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  4;
		   ind_loop[2]=  1- param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1]+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);
		}
	       }
	}
   





  

  //delete [] iptdrodm;
 delete [] param_int;
   

  //RELEASESHAREDN( stock  , stk );
  //RELEASESHAREDN( drodmArray,drodm);
  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
  
} 
