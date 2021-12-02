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
  PyObject* constk;
  E_Int taille_tabs;
 
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOOOll", &zones , &stock, &drodmstock, &constk, &work, &nstep, &taille_tabs)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOOOii", &zones , &stock, &drodmstock, &constk, &work, &nstep, &taille_tabs)) return NULL; 
#endif

  vector<PyArrayObject*> hook;

  E_Int nidom = PyList_Size(zones);

  E_Int**   param_int = new E_Int*[nidom];
  E_Float** param_real= new E_Float*[nidom];
  E_Float** iptro_p1 = new E_Float*[nidom];
  E_Float** iptro    = new E_Float*[nidom];
  E_Float** iptro_m1    = new E_Float*[nidom];

  
  /// Recuperation du tableau de stockage des valeurs
  FldArrayF* stk;
  K_NUMPY::getFromNumpyArray(stock, stk, true); E_Float* iptstk = stk->begin();

  /// Recuperation du tableau de stockage des flux
  FldArrayF* drodmstk;
  K_NUMPY::getFromNumpyArray(drodmstock,drodmstk, true); E_Float* iptdrodmstk = drodmstk->begin();

  /// Recuperation du tableau de stockage des flux pour conservativite
  FldArrayF* cstk;
  K_NUMPY::getFromNumpyArray(constk, cstk, true); E_Float* iptcstk = cstk->begin();

  /// Recuperation du tableau des flux (drodm)
  PyObject* drodmArray = PyDict_GetItemString(work,"rhs"); FldArrayF* drodm;
  K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();
   
  /// Tableau de travail coe   ( dt/vol et diags LU)
  PyObject* coeArray = PyDict_GetItemString(work,"coe"); FldArrayF* coe;
  K_NUMPY::getFromNumpyArray(coeArray, coe, true); E_Float* iptcoe = coe->begin();

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

       o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t2  = K_PYTREE::getNodeFromName1( o        , "Density_M1");
       iptro_m1[nd]     = K_PYTREE::getValueAF(t1, hook);
      }  

  E_Int imax=0;
  E_Int c=0;
  E_Int shift;
  E_Int cycle;
  E_Int a;
  a=0;
  E_Int b;
  b=0;
  E_Int shift_zone[nidom];
  E_Int shift_coe[nidom];

     for (E_Int nd = 0; nd < nidom; nd++)
       {
	 shift_zone[nd]=a;
	 a=a+param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
	 //cout << "ndimdx= " << param_int[nd][ NDIMDX ] << endl;
       }
      for (E_Int nd = 0; nd < nidom; nd++)
       {
	 shift_coe[nd]=b;
	 b=b+param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];	 
       }



   if (param_int[0][RK]==2 and param_int[0][EXPLOC]==1) // Mon schéma
      {      
      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  if (nstep%(param_int[0][NSSITER]/param_int[nd][LEVEL])==1) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	      
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// stockage des flux
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		  
		   /// stockage des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   /// Interpolation
		   interpolation4_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);

		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// stockage des flux
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;

		   /// stockage des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  2;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   /// Interpolation
		   interpolation4_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);

		}
	     
	    }    



	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==(param_int[0][NSSITER]/param_int[nd][LEVEL])/2 and (param_int[0][NSSITER]/param_int[nd][LEVEL]) !=  2) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     

	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   /// Conservativite : Stockage du flux
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV];
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   		   
		   for (E_Int n = 0; n < nd+1 ; n++)
		     //for (E_Int n = 0; n < 0 ; n++)
		     {
		       c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		    shift=c;
		    //cout <<"shift_coe= "<< shift_coe[nd] << endl;

		    c=nd+1;
		    //c = nd+2;
		    //c=0;
		   //cout <<"c=  "<< c << endl;
		   imax=1;



		   conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0; 

		   //cout << "coucou conservativite 3" << endl;
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		   /// Conservativite : Modif du bilan des flux sur la premiere maille du domaine
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  1;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   	      
		   for (E_Int n = 0; n < nd-1 ; n++)
		   //for (E_Int n = 0; n < 0 ; n++)
		      {
		              c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		   
		   shift=c;
		   //cout <<"c(cons2)=  "<< c << endl;
		   //c=nd-1;
		   c=nd-1;
		   //c=0;
		   imax=param_int[c][IJKV];
		   //cout <<"param(IJKV)(cons2)=  "<< imax << endl;

		   //cout << "coucou conservativite 3" << endl;

		   conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0;

		   //cout << "coucou conservativite 3" << endl;
		   
		}
	    }



	  if (((nstep+1)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// Recuperation des flux
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		  //  ind_loop[0]=  param_int[nd][IJKV]-2;
		  //  ind_loop[1]=  param_int[nd][IJKV];
		  //  ind_loop[2]=  1;
		  //  ind_loop[3]=  param_int[nd][IJKV+1];
		  //  ind_loop[4]=  1;
		  //  ind_loop[5]=  param_int[nd][IJKV+2];
       
		  //  copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);

		  
		   /// Recuperation des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// Recuperation des flux
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		  //  ind_loop[0]=  1;
		  //  ind_loop[1]=  3;
		  //  ind_loop[2]=  1;
		  //  ind_loop[3]=  param_int[nd][IJKV+1];
		  //  ind_loop[4]=  1;
		  //  ind_loop[5]=  param_int[nd][IJKV+2];
       
		  //  copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);

		   /// Recuperation des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  2;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);
		}
	       }
    

	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and (param_int[0][NSSITER]/param_int[nd][LEVEL]) !=  2) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   /// Conservativite : Stockage du flux
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV];
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   		   
		   for (E_Int n = 0; n < nd+1 ; n++)
		   //for (E_Int n = 0; n < 0 ; n++)
		      {
		       c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		    shift=c;
		    //cout <<"c=  "<< c << endl;
		    c=nd+1;
		    //c=0;
		   imax=1;
		   //conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0;  
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		   /// Conservativite : Modif du bilan des flux sur la premiere maille du domaine
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  1;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   	      
		   for (E_Int n = 0; n < nd-1 ; n++)
		   //for (E_Int n = 0; n < 0 ; n++)
		      {
		              c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		   
		   shift=c;
		   //cout <<"c(cons2)= "<< c << endl;
		   c=nd-1;
		   //c=0;
		   imax=param_int[c][IJKV];
		   //cout <<"param(IJKV)(cons2)= "<< imax << endl;
		 
		   //conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0;
		   
		}
	    }
	  
	  	  
	}
  
      }








   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////TANG & WARNECKE///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







   if (param_int[0][RK]==2 and param_int[0][EXPLOC]==5) //Schema de Tang & Warnecke
      {      
      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  if (nstep%(param_int[0][NSSITER]/param_int[nd][LEVEL])==1) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	      
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
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
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   /// Interpolation
		   interpolation2_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);

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
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   /// Interpolation
		   interpolation2_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);

		}
	     
	    }    



	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==(param_int[0][NSSITER]/param_int[nd][LEVEL])/2 and (param_int[0][NSSITER]/param_int[nd][LEVEL]) !=  2) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   /// Conservativite : Stockage du flux
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV];
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   		   
		   for (E_Int n = 0; n < nd+1 ; n++)
		     {
		       c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		    shift=c;
		    //cout <<"shift_coe= "<< shift_coe[nd] << endl;

		   c=nd+1;
		   //cout <<"c=  "<< c << endl;
		   imax=1;

		   //cout << "coucou conservativite 3" << endl;

		   conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0; 

		   //cout << "coucou conservativite 3" << endl;
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		   /// Conservativite : Modif du bilan des flux sur la premiere maille du domaine
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  1;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   	      
		   for (E_Int n = 0; n < nd-1 ; n++)
		      {
		              c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		   
		   shift=c;
		   //cout <<"c(cons2)=  "<< c << endl;
		   c=nd-1;
		   imax=param_int[c][IJKV];
		   //cout <<"param(IJKV)(cons2)=  "<< imax << endl;

		   //cout << "coucou conservativite 3" << endl;

		   conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0;

		   //cout << "coucou conservativite 3" << endl;
		   
		}
	    }



	  if (((nstep+1)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
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
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

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
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);
		}
	       }
    

	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and (param_int[0][NSSITER]/param_int[nd][LEVEL]) !=  2) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   /// Conservativite : Stockage du flux
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV];
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   		   
		   for (E_Int n = 0; n < nd+1 ; n++)
		      {
		       c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		    shift=c;
		    //cout <<"c=  "<< c << endl;
		   c=nd+1; 
		   imax=1;
		   //conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0;  
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		   /// Conservativite : Modif du bilan des flux sur la premiere maille du domaine
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  1;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   	      
		   for (E_Int n = 0; n < nd-1 ; n++)
		      {
		              c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		      }
		   
		   shift=c;
		   //cout <<"c(cons2)= "<< c << endl;
		   c=nd-1;
		   imax=param_int[c][IJKV];
		   //cout <<"param(IJKV)(cons2)= "<< imax << endl;
		 
		   //conservativite32_(idir,param_int[nd],param_real[nd],param_int[c],ind_loop,imax,iptro[nd],iptdrodm + shift_zone[nd],iptdrodm + shift,iptcstk,iptcoe+shift_coe[nd],nstep,nd);
		   c=0;
		   
		}
	    }
	  
	  	  
	}
  
      }







   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////CONSTANTINESCU RK2///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







      if (param_int[0][RK]==2 and param_int[0][EXPLOC]==2) // Schéma de Constantinescu
      {      
      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  if (nstep%(param_int[0][NSSITER]/param_int[nd][LEVEL])==1 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {

	       if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// stockage des flux (drodm)
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  param_int[nd][IJKV]-4;
		   ind_loop[1]=  param_int[nd][IJKV]-4;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

		   copyflux_( idir, param_int[nd], ind_loop,iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);
		   		  		  
		   /// stockage des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-3;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

       		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptstk,ind,nd);		   
	   
		  		   
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{

		  //// stockage des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  5;
		   ind_loop[1]=  5;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);

		   /// stockage des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  4;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptstk,ind,nd);

   		
		}	     
	    }    

	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==(param_int[0][NSSITER]/param_int[nd][LEVEL])/2 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;

		   /// Recuperation des valeurs (yn)
		   ind_loop[0]=  param_int[nd][IJKV]-3;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

       		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptstk,ind,nd);

		   /// Stockage du flux à ajouter à la derniere etape
		   ind=1;
		   ind_loop[0]=  param_int[nd][IJKV]-3;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   constantinescu_(iptro[nd],ind_loop,iptcstk,iptdrodm+shift_zone[nd],iptcoe+shift_coe[nd],param_int[nd],param_real[nd],nstep,nd);


		}
	     
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
	
		{

		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   /// Recuperation des valeurs (yn)
		   ind_loop[0]=  1;
		   ind_loop[1]=  4;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptstk,ind,nd); 

		   /// Stockage du flux à ajouter à la derniere etape
		   ind=1;
		   ind_loop[0]=  1;
		   ind_loop[1]=  4;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   constantinescu_(iptro[nd],ind_loop,iptcstk,iptdrodm +shift_zone[nd],iptcoe+shift_coe[nd],param_int[nd],param_real[nd],nstep,nd);

		   
		}
	    }
	  
       

	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==(param_int[0][NSSITER]/param_int[nd][LEVEL])/2+1 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// Recuperation des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  param_int[nd][IJKV]-4;
		   ind_loop[1]=  param_int[nd][IJKV]-4;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);   		   
	  
	    
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// Recuperation des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  5;
		   ind_loop[1]=  5;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);
	   
		   
		}
	    }


	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{


		  //cout << nd << endl;

		   /// Derniere etape schema de Constantinescu
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV]-3;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   		   
		   // for (E_Int n = 0; n < param_int[nd][NUMD] ; n++)
		   //   {
		   //    c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		   //   }
		   // shift=c;
		   //c=param_int[nd][NUMD]; 
		   //imax=1;
		   //cout << "shift_zone=  " << shift_zone[nd] << endl;
		   
		   constantinescu_(iptro[nd],ind_loop,iptcstk,iptdrodm +shift_zone[nd],iptcoe+shift_coe[nd],param_int[nd],param_real[nd],nstep,nd);
		   //c=0;  
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  

		  //cout << nd << endl;

		   /// Derniere etape schema de Constantinescu
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  4;
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   	      
		   //for (E_Int n = 0; n < param_int[nd][NUMG] ; n++)
		   //   {
		   //           c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		   //   }
		   
		   //shift=c;
		   //cout <<"c(cons2)= "<< c << endl;
		   //c=param_int[nd][NUMG];
		   //imax=param_int[c][IJKV];
		   //cout <<"param(IJKV)(cons2)= "<< imax << endl;

		   //cout << "shift_zone=  " << shift_zone[nd] << endl;
		 
		   constantinescu_(iptro[nd],ind_loop,iptcstk,iptdrodm +shift_zone[nd],iptcoe+shift_coe[nd],param_int[nd],param_real[nd],nstep,nd);
		   //c=0;
		   
		}
	    }
	  
	  	  
	}
      }









   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////CONSTANTINESCU RK3//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


     



      if (param_int[0][RK]==2 and param_int[0][EXPLOC]==3) // Schéma de Constantinescu base sur rk3
      {      
      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  if (nstep%(param_int[0][NSSITER]/param_int[nd][LEVEL])==2 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	       if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //cout <<ipt_param_int[nd][IJKV]-2 <<endl;
		  //cout <<ipt_param_int[nd][IJKV]   <<endl;
		  //cout <<ipt_param_int[nd][IJKV+1] <<endl;
		  //cout <<"nzone = "<<nd<<"stockage flux cel adj slow buffer et valeurs"<<endl;
		  
		  //// stockage du drodm de la maille adjacente au slow buffer
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  param_int[nd][IJKV]-6;
		   ind_loop[1]=  param_int[nd][IJKV]-6;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   //cout << "shift_zone[nd]= " <<shift_zone[nd]<< endl; 
		   copyflux_( idir, param_int[nd], ind_loop,iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd); 

		   /// Stockage des valeurs y^n 
		   ind_loop[0]=  param_int[nd][IJKV]-5;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

       		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptcstk,ind,nd);		  		  
   
	   
		  		   
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //cout << "stockage,interpolation,conservativite"<<endl;
		  //// stockage des flux
		  //cout <<"nzone = "<<nd<<"stockage flux cel adj slow buffer et valeurs"<<endl;


		  //// stockage du drodm de la maille adjacente au slow buffer
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  7;
		   ind_loop[1]=  7;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);

		   /// Stockage des valeurs y^n 
		   ind_loop[0]=  1;
		   ind_loop[1]=  6;
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

      		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptcstk,ind,nd); 

   		
		}	     
	    }    

	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==(param_int[0][NSSITER]/param_int[nd][LEVEL])/2 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   //cout <<"coucou"<<endl;

		   /// stockage des valeurs y^3
		   ind_loop[0]=  param_int[nd][IJKV]-5;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

       		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptstk,ind,nd);		   		   
		   
		   ind=2;
		   /// Recuperation des valeurs y^n ds rop
		   ind_loop[0]=  param_int[nd][IJKV]-5;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

       		   copy_(idir, param_int[nd], ind_loop, iptro[nd], iptcstk, ind, nd);


		}
	     
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
	
		{
		   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   //cout <<"nzone = "<<nd<<"recup valeurs et stockage flux derniere etape"<<endl;	
		   //cout <<"coucou"<<endl;
		   /// stockage des valeurs y^3
		   ind_loop[0]=  1;
		   ind_loop[1]=  6;
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro[nd], iptstk,ind,nd);


		   ind=2;
         	   /// Recuperation des valeurs y^n ds rop
		   ind_loop[0]=  1;
		   ind_loop[1]=  6;
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
       		   copy_(idir, param_int[nd], ind_loop, iptro[nd], iptcstk, ind, nd);	   
		}
	    }
	  
       

	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==(param_int[0][NSSITER]/param_int[nd][LEVEL])-1 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //cout <<ipt_param_int[nd][IJKV]-2 <<endl;
		  //cout <<ipt_param_int[nd][IJKV]   <<endl;
		  //cout <<ipt_param_int[nd][IJKV+1] <<endl;
		  //cout <<ipt_param_int[nd][IJKV+2] <<endl;
		  //cout <<"nzone = "<<nd<<"recuperation"<<endl;
		  //cout <<"nzone = "<<nd<<"recup flux cel adj slow buffer"<<endl;
			  
		  //// Recuperation du drodm maille adjacente slow buffer
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  param_int[nd][IJKV]-6;
		   ind_loop[1]=  param_int[nd][IJKV]-6;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);   		   
	  
	    
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{

		  // cout <<"nzone = "<<nd<<"recup flux cel adj slow buffer"<<endl;

		  //// Recuperation drodm maillae adjacente slow buffer
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  7;
		   ind_loop[1]=  7;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);
	   
		   
		}
	    }


	  if (((nstep)%(param_int[0][NSSITER]/param_int[nd][LEVEL]))==0 and param_int[nd][LEVEL]<(param_int[0][NSSITER]/param_int[nd][RK])) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		   /// Derniere etape schema de Constantinescu
		  //cout <<"nzone = "<<nd<<"constantinescu derniere etape"<<endl;	
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV]-5;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   		   
		   constantinescurk32_(idir,param_int[nd],param_real[nd],ind_loop,iptro[nd],iptstk,nd);
  
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		   /// Derniere etape schema de Constantinescu
		  //cout <<"nzone = "<<nd<<"constantinescu derniere etape"<<endl;	
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  1;
		   ind_loop[1]=  6;
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
		   	      
		   //for (E_Int n = 0; n < param_int[nd][NUMG] ; n++)
		   //   {
		   //           c=c+param_int[n][ NDIMDX ]*param_int[n][ NEQ ];	 
		   //   }
		   
		   //shift=c;
		   //cout <<"c(cons2)= "<< c << endl;
		   //c=param_int[nd][NUMG];
		   //imax=param_int[c][IJKV];
		   //cout <<"param(IJKV)(cons2)= "<< imax << endl;


		   constantinescurk32_(idir,param_int[nd],param_real[nd],ind_loop,iptro[nd],iptstk,nd);
		 

		   //c=0;
		   
		}
	    }
	  
	  	  
	}
      }





     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////RK3 LOCAL/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



   E_Int taille=2000000/nidom;
   //cout << nidom << endl;
   E_Int k=1;

     
   if(param_int[0][RK]==3 and param_int[0][EXPLOC]==2)
     {
       for (E_Int nd = 0; nd < nidom; nd++)
	 {
	   cycle=param_int[0][NSSITER]/param_int[nd][LEVEL];
	  
	   if (cycle==4) //////////////////////////////////// Zones de plus petit niveau en temps ////////////////////////////////////////////////// 
	     {

	       //cout << nd << endl;

	       if (nstep%cycle==cycle/2-1 and (nstep/cycle)%2==0) /// Stockage pour la zone adjacente de + petit niveau en tps
		 {
		   if(param_int[nd][LEVELG]==0) /// La zone de + gd niveau en temps est sur le bord gauche du domaine : elle fournit des infos à droite
		     {

		       //// stockage des flux en vue de l'interpolation pour la zone de + petit niveau en tps
		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;E_Int pos=0; 
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
		       //cout << "coucou" << endl;
		       copyfluxrk3local_( idir, param_int[nd], ind_loop,iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		       		  
		       /// stockage de yn en position 0 dans le tableau de stockage
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
		       //cout << "coucou" << endl;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  


		     }
		   else if(param_int[nd][LEVELD]==0) /// La zone de + gd niveau en temps est sur le bord droit du domaine : elle fournit des infos à gauche
		     {
		       //// stockage des flux en vue de l'interpolation pour la zone de + petit niveau en tps
		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;E_Int pos=0;
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		       		  
		       /// stockage de yn en position 0 dans le tableau de stockage
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);
		   		   
		     }	  

	           else  /// La zone de + gd niveau en temps est entourée par 2 zones : elle fournit des infos des deux côtés

		     {

		       //// stockage des flux en vue de l'interpolation pour la zone de + petit niveau en tps
		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1; E_Int pos=0; 
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       //cout << "shift_zone= " << shift_zone[1]<< endl;

		       copyfluxrk3local_( idir, param_int[nd], ind_loop,iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		      		  
		       /// stockage de yn en position 0 dans le tableau de stockage
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille); 

		       //// stockage des flux en vue de l'interpolation pour la zone de + petit niveau en tps
		       pos=1;
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];

		       //cout << "shift_zone=  " <<shift_zone[1] << endl;
       
		       copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		       		  
		       /// stockage de yn en position 1 dans le tableau de stockage
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
		       

       
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

		     }




		 }    

	       if (nstep%cycle==cycle/2 ) /// Recuperation des valeurs stockées en position 0 (yn ou yinterpolé suivant la parité du cycle) et stockage de y2 ou y6
		 {
	     
		   if(param_int[nd][LEVELG]==0) /// savoir quelles colonnes on prend
		     {
	
		       E_Int idir = 2; E_Int ind_loop[6]; 		  

		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
		       
		       if ((nstep/cycle)%2==1)
			 {
			   /// Stockage de y2 ou y6 (suivant la parité du cycle) en position 1 dans le tableau de stockage
			   E_Int ind=1;E_Int pos=1;
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			   /// Recuperation en position 0 dans le tableau de stockage de yinterpolé
			   ind=2;pos=0;       
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  
			 }

		        if ((nstep/cycle)%2==0)
		        {
			  //cout << "coucou" << endl;
			  E_Int ind=2;E_Int pos=0;
		         copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);		       
			 }



		     }
		   else if(param_int[nd][LEVELD]==0) /// savoir quelles colonnes on prend
		     {

		       E_Int idir = 2; E_Int ind_loop[6]; 
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];      
		       
		       if ((nstep/cycle)%2==1)
			 {
			   /// Stockage de y6 en position 1 dans le tableau de stockage
			   E_Int ind=1;E_Int pos=1;
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			   /// Recuperation en position 0 dans le tableau de stockage de yinterpolé
			   ind=2;pos=0;       
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  
			 }

		       if ((nstep/cycle)%2==0)
		        {
			  //cout << "coucou" <<" "<<pos<< endl;
		          E_Int ind=2;E_Int pos=0;
		          copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);
		        }

		   		   
		     }	 
		   else
		     {

		       E_Int idir = 2; E_Int ind_loop[6];E_Int ind=1;E_Int pos=2; 		  

		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];

		       if ((nstep/cycle)%2==1)
			 {
			   /// Stockage des 2 dernieres colonnes de  y6 en position 2 dans le tableau de stockage
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			   /// Recuperation en position 0 dans le tableau de stockage des 2 dernieres colonnes de yinterpolé
			   ind=2;pos=0;       
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			   pos=3;ind=1; 
			   ind_loop[0]=  1;
			   ind_loop[1]=  2;

			   /// Stockage des 2 premieres colonnes de y6 en position 3 dans le tableau de stockage
		           copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			   /// Recuperation en position 1 dans le tableau de stockage de 2 premieres colonnes  de yinterpolé
			   ind=2; pos=1;       
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  
			  
			 }
			   
		        if ((nstep/cycle)%2==0)
			  {

			  ind_loop[0]=  param_int[nd][IJKV]-1;
			  ind_loop[1]=  param_int[nd][IJKV];

			  E_Int ind=2; E_Int pos=0;
		          copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

			  ind_loop[0]=  1;
			  ind_loop[1]=  2;

			  ind=2; pos=1;
		          copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);		       
			  }

		     }



		 }



	       if (nstep%cycle==0) /// Interpolation pour la zone adjacente de pas de temps plus petit et switch pointeurs
		 {

		   // Switch pointeurs
		   switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;
		   
	   
		   if ((nstep/cycle)%2==1 )
		     {
		  
		       if(param_int[nd][LEVELG]==0) /// savoir quelles colonnes on prend
			 {
			   //cout << "coucou1" << endl;
			   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;E_Int pos=0;		  
			   /// Recuperation en position 0 dans le tableau de stockage
			   ind_loop[0]=  param_int[nd][IJKV]-1;
			   ind_loop[1]=  param_int[nd][IJKV];
			   ind_loop[2]=  1;
			   ind_loop[3]=  param_int[nd][IJKV+1];
			   ind_loop[4]=  1;
			   ind_loop[5]=  param_int[nd][IJKV+2];
       
			   interprk3local_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptdrodmstk+nd*taille, nd, pos,taille);  


			 }
		       else if(param_int[nd][LEVELD]==0) /// savoir quelles colonnes on prend
			 {
			   //cout << "coucou2" << endl;
			   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2; E_Int pos=0;
			   /// Recuperation de yn en position 0 dans le tableau de stockage
			   //cout << "coucou" <<" "<<pos<< endl;


			   ind_loop[0]=  1;
			   ind_loop[1]=  2;
			   ind_loop[2]=  1;
			   ind_loop[3]=  param_int[nd][IJKV+1];
			   ind_loop[4]=  1;
			   ind_loop[5]=  param_int[nd][IJKV+2];
       
			   interprk3local_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptdrodmstk+nd*taille, nd, pos,taille);  
		       		   		   
			 }
		       else
			 {
			   E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;E_Int pos=0;		  
			   /// Recuperation en position 0 dans le tableau de stockage
			   ind_loop[0]=  param_int[nd][IJKV]-1;
			   ind_loop[1]=  param_int[nd][IJKV];
			   ind_loop[2]=  1;
			   ind_loop[3]=  param_int[nd][IJKV+1];
			   ind_loop[4]=  1;
			   ind_loop[5]=  param_int[nd][IJKV+2];
       
			   interprk3local_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptdrodmstk+nd*taille, nd, pos,taille);

			   pos=1;
			   /// Recuperation en position 1 ds le tableau de stockage
			   ind_loop[0]=  1;
			   ind_loop[1]=  2;
			   ind_loop[2]=  1;
			   ind_loop[3]=  param_int[nd][IJKV+1];
			   ind_loop[4]=  1;
			   ind_loop[5]=  param_int[nd][IJKV+2];  

			   interprk3local_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptdrodmstk+nd*taille, nd, pos,taille);  
			 }
		     }
		 }
	     
	     }

	   else if (cycle == param_int[nd][NSSITER]) /// Zones de + gd pas de temps 

	     {

	       //cout << nd << endl;


	       if (param_int[nd][LEVELG]==0) /// La zone est situee sur le bord gauche du domaine et transmet des infos à droite
		 {
		   E_Int idir = 2; E_Int ind_loop[6];
		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];		   
		    
		   if (nstep%cycle == 1)
		     {

		       /// Stockage de yn en position 0 (4 dernieres colonnes)
		       ind_loop[0]=  param_int[nd][IJKV]-3;
		       ind_loop[1]=  param_int[nd][IJKV];
		       E_Int ind=1;E_Int pos=0;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  


		       /// Stockage de y3 (4 dernieres colonnes) en position 2 car il y a deja les 4 dernieres colonnes de yn
		       pos=2;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille); 


		       /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn (sur les 4 dernieres colonnes)
		       pos=2;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);

		       /// Stockage des 3 dernieres colonnes du drodm
		       ind_loop[0]=  param_int[nd][IJKV]-2;
		       ind_loop[1]=  param_int[nd][IJKV]+1;
		       ind=1;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);


		     }


		   if (nstep%cycle == cycle/2-1)
		     {
		       /// Recuperation de y3 (4 dernieres colonnes)
		       ind_loop[0]=  param_int[nd][IJKV]-3;
		       ind_loop[1]=  param_int[nd][IJKV];
		       E_Int ind=2;E_Int pos=2;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);

		       /// Recuperation des 3 dernieres colonnes du drodm
		       ind_loop[0]=  param_int[nd][IJKV]-2;
		       ind_loop[1]=  param_int[nd][IJKV]+1;
		       ind=2;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);
		     }


		   if (nstep%cycle == cycle/2)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];       
      		      
		       /// Stockage de y4 en position 4 car il y a deja les 4 dernieres colonnes de yn et les 4 dernieres colonnes de y1
		       E_Int ind=1;E_Int pos=4;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);
		     }

		   if (nstep%cycle == cycle/2 + cycle/4 - 1)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       E_Int ind=1;E_Int pos=2;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);
		     }
		   if (nstep%cycle == cycle/2 + cycle/4) 
		     {
		       /// Interpolation de y6
		       E_Int ind=1;E_Int pos=2;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro[nd], nd, pos,taille,nstep);
		     }

		   if (nstep%cycle == 0) 
		     {

		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;
	      
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
		    
		   if (nstep%cycle == 1)
		     {
		       /// Stockage de yn en position 0 (4 premieres colonnes)
		       ind_loop[0]= 1;
		       ind_loop[1]= 4; 
		       E_Int ind=1;E_Int pos=0;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  


		       /// Stockage de y3 en position 2 car il y a deja les 4 premieres colonnes de yn
		       pos=2;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille); 


		       /// Interpolation de y1
		       pos=2;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);


		       /// Stockage des 3 premieres colonnes du drodm
		       ind_loop[0]=  0;
		       ind_loop[1]=  3;
		       ind=1;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);
 
		     }

		   if (nstep%cycle == cycle/2-1)
		     {
		       /// Recuperation de y3 (4 premieres colonnes)
		       ind_loop[0]=  1;
		       ind_loop[1]=  4;
		       E_Int ind=2;E_Int pos=2;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);

		       /// Recuperation des 3 premieres colonnes du drodm
		       ind_loop[0]=  0;
		       ind_loop[1]=  3;
		       ind=2;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);
		     }

		   if (nstep%cycle == cycle/2)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       ind_loop[0]=  1;
		       ind_loop[1]=  2;

		       /// Stockage de y4 en position 4
		       E_Int ind=1;E_Int pos=4;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);
		     }

		   if (nstep%cycle == cycle/2 + cycle/4 - 1)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       /// Interpolation de y5
		       E_Int ind=1;E_Int pos=2;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);
		     }
		   if (nstep%cycle == cycle/2 + cycle/4) 
		     {
		       /// Interpolation de y6
		       E_Int ind=1;E_Int pos=2;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro[nd], nd, pos,taille,nstep);
		     }

		   if (nstep%cycle == 0) 
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;
	      
		     }



		 }
	     
	       
	     }
	   else
	     {
	       if (param_int[nd][LEVELG] < param_int[nd][LEVEL] and param_int[nd][LEVEL] < param_int[nd][LEVELD])

		 //cout << "coucou" << endl;

		 {

		       E_Int idir = 2; E_Int ind_loop[6];
		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
 
		       //cout << "coucou" << endl;

		   if (nstep%cycle == 1)
		     {

		       ind_loop[0]=  param_int[nd][IJKV]-3;
		       ind_loop[1]=  param_int[nd][IJKV];

		       /// Stockage des 4 dernieres colonnes de yn en position 0 
		       E_Int ind = 1; E_Int pos=0;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

		       ind_loop[0]=  param_int[nd][IJKV]-3;
		       ind_loop[1]=  param_int[nd][IJKV];

		       /// Stockage des 4 dernieres colonnes de y3 en position 6 
		       ind = 1; pos=3;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);

		       /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn (sur les 4 dernieres colonnes)
		       pos=3;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);

		       /// Stockage des 3 dernieres colonnes du drodm en position 0
		       ind_loop[0]=  param_int[nd][IJKV]-2;
		       ind_loop[1]=  param_int[nd][IJKV]+1;
		       ind=1;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		       if ((nstep/cycle)%2==0)
			 {

		           ind_loop[0]=  1;
			   ind_loop[1]=  2;

			   /// Stockage des 2 premieres colonnes de yn en position 4 
			   E_Int ind = 1; E_Int pos=2;
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

			   /// Stockage des 2 premieres colonnes du drodm en position 4
			   ind=1;pos=2;
			   copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);
			 }

		     }

		   if (nstep%cycle == cycle/2 - 1)
		     {

		       /// Recuperation de y3 (4 dernieres colonnes)
		       ind_loop[0]=  param_int[nd][IJKV]-3;
		       ind_loop[1]=  param_int[nd][IJKV];
		       E_Int ind=2;E_Int pos=3;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);

		       /// Recuperation des 3 dernieres colonnes du drodm
		       ind_loop[0]=  param_int[nd][IJKV]-2;
		       ind_loop[1]=  param_int[nd][IJKV]+1;
		       ind=2;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		     }

		   if (nstep%cycle == cycle/2)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];
     		      
		       /// Stockage des 2 dernieres colonnes de y4 en position 10 
		       E_Int ind=1;E_Int pos=5;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille); 

		       if ((nstep/cycle)%2==0)
		        {

			  //// Recuperation et transformation de f(yn) + coeff*f(y1) pour obtenir alpha*f(yn) + beta*f(y1)
			  ind_loop[0]=  1;
			  ind_loop[1]=  2;
			  ind=2; pos=2;
			  copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		        }


		     }

		   if (nstep%cycle == cycle/2 + cycle/4-1)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       ind_loop[0]=  param_int[nd][IJKV]-1;
		       ind_loop[1]=  param_int[nd][IJKV];

		       /// Interpolation de y5
		       E_Int ind=1;E_Int pos=3;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);

		     }


		   if (nstep%cycle == cycle/2 + cycle/4)
		     {

		       /// Interpolation de y6
		       E_Int ind=1;E_Int pos=3;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro[nd], nd, pos,taille,nstep);

		     }

		   if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1 )

		        {

			  //// Stockage des 2 premieres colonnes de y6 en position 12
			  ind_loop[0]=  1;
			  ind_loop[1]=  2;
			  E_Int ind=1; E_Int pos=6;
			  copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			  //// Recuperation de yinterpolé en position 4
			  ind=2; pos=2;
			  copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

		        }


		   if (nstep%cycle == 0)
		     {

		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;


		       if ((nstep/cycle)%2==1 )
			 {

			   ind_loop[0]=  1;
			   ind_loop[1]=  2;

			   E_Int pos=2;
			   interprk3local_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptdrodmstk+nd*taille, nd, pos,taille);  

			 }


		     }



	       
		     }



	       if (param_int[nd][LEVELG] > param_int[nd][LEVEL] and param_int[nd][LEVEL] > param_int[nd][LEVELD])

		 //cout << nd << endl;

		 {

		       E_Int idir = 2; E_Int ind_loop[6];
		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];
 
		       //cout << "coucou" << endl;

		   if (nstep%cycle == 1)
		     {

		       ind_loop[0]=  1;
		       ind_loop[1]=  4;

		       /// Stockage des 4 premieres colonnes de yn en position 0 
		       E_Int ind = 1; E_Int pos=0;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

		       ind_loop[0]=  1;
		       ind_loop[1]=  4;

		       /// Stockage des 4 premieres colonnes de y3 en position 6 
		       ind = 1; pos=3;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);


		       /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn (sur les 4 premieres colonnes)
		       pos=3;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);

		       /// Stockage des 3 premieres colonnes du drodm en position 0
		       ind_loop[0]=  0;
		       ind_loop[1]=  3;
		       ind=1;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);


		       if ((nstep/cycle)%2==0)
			 {

			   ind_loop[0]=  param_int[nd][IJKV]-1;
			   ind_loop[1]=  param_int[nd][IJKV];

			   /// Stockage des 2 dernieres colonnes de yn en position 4 
			   E_Int ind = 1; E_Int pos=2;
			   copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

			   /// Stockage des 2 dernieres colonnes du drodm en position 4
			   pos=2;ind=1;
			   ind_loop[0]= param_int[nd][IJKV]-1;
			   ind_loop[1]= param_int[nd][IJKV];
			   copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

			 }


		     }

		   if (nstep%cycle == cycle/2 - 1)
		     {
		       /// Recuperation de y3 (4 premieres colonnes)
		       ind_loop[0]=  1;
		       ind_loop[1]=  4;
		       E_Int ind=2;E_Int pos=3;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille);

		       /// Recuperation des 3 premieres colonnes du drodm
		       ind_loop[0]=  0;
		       ind_loop[1]=  3;
		       ind=2;pos=0;
		       copyfluxrk3local2_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);
  		    

		     }

		   if (nstep%cycle == cycle/2)

		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       ind_loop[0]=  1;
		       ind_loop[1]=  2;
     		      
		       /// Stockage des 2 premieres colonnes de y4 en position 10 
		       E_Int ind=1;E_Int pos=5;
		       copyrk3local_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk + nd*taille,ind,pos,taille); 

		       if ((nstep/cycle)%2==0)
		        {

			  //// Recuperation et transformation de f(yn) + coeff*f(y1) pour obtenir alpha*f(yn) + beta*f(y1)
			  ind_loop[0]= param_int[nd][IJKV]-1;
			  ind_loop[1]= param_int[nd][IJKV];
			  ind=2; pos=2;
			  copyfluxrk3local_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk+nd*taille,ind,pos,taille);

		        }

		     }

		   if (nstep%cycle == cycle/2 + cycle/4-1)
		     {
		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;

		       ind_loop[0]=  1;
		       ind_loop[1]=  2;

		       /// Interpolation de y5
		       E_Int ind=1;E_Int pos=3;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro_p1[nd], nd, pos,taille,nstep);

		     }


		   if (nstep%cycle == cycle/2 + cycle/4)
		     {

		       /// Interpolation de y6
		       E_Int ind=1;E_Int pos=3;
		       interprk3local2_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptro[nd], nd, pos,taille,nstep);

		     }


		   if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1 )

		    {

			  //// Stockage des 2 dernieres colonnes de y6 en position 12
			  ind_loop[0]=  param_int[nd][IJKV]-1;
			  ind_loop[1]=  param_int[nd][IJKV] ;
			  E_Int ind=1; E_Int pos=6;
			  copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);  

			  //// Recuperation de yinterpolé en position 4
			  ind=2; pos=2;
			  copyrk3local_( idir, param_int[nd], ind_loop, iptro[nd], iptstk + nd*taille,ind,pos,taille);

		    }



		   if (nstep%cycle == 0)
		     {

		       /// Switch pointeurs
		       switchvectors_(iptro[nd],iptro_p1[nd],param_int[nd]) ;


		       if ((nstep/cycle)%2==1 )
			 {

			   ind_loop[0]= param_int[nd][IJKV]-1 ;
			   ind_loop[1]= param_int[nd][IJKV] ;

			   E_Int pos=2;
			   interprk3local_(idir, param_int[nd], param_real[nd], iptcoe+shift_coe[nd], ind_loop, iptstk + nd*taille, iptdrodmstk+nd*taille, nd, pos,taille);  

			 }


		     }




		       


		     }





		 }


	     }

       k=2;
	   
     }
    




   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////RK2 LOCAL TEST//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




     if (param_int[0][RK]==2 and param_int[0][EXPLOC]==4) // RK2 local test
      {      
      for (E_Int nd = 0; nd < nidom; nd++)
	{
	  cycle=param_int[0][NSSITER]/param_int[nd][LEVEL];
	  
	  if (nstep%cycle==1 and param_int[nd][LEVEL] < cycle) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {

	       if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// stockage des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  param_int[nd][IJKV]-4;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop,iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);
		   
		   		  		  
		   /// stockage des valeurs
		   ind_loop[0]=  param_int[nd][IJKV]-5;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];


		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,nd);

		   
		   /// Interpolation
		   interpolation2_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]);

		  
	   
		  		   
		}
	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{

		  //// stockage des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;
		   ind_loop[0]=  1;
		   ind_loop[1]=  5;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm + shift_zone[nd], iptdrodmstk,ind,nd);


		   /// stockage des valeurs
		   ind_loop[0]=  1;
		   ind_loop[1]=  6;
		   ind_loop[2]=  1;//-param_int[nd][NIJK+3];
		   ind_loop[3]=  param_int[nd][IJKV+1];//+param_int[nd][NIJK+3];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

		   E_Int pos=1;       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);

		   /// Interpolation
		   interpolation2_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd]); 


	   
		
		}	     
	    }    

	if (nstep==1 and param_int[nd][LEVEL] == cycle) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {

	      //cout << "coucou" << endl;
	       if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1; E_Int pos=6;
		   ind_loop[0]=  1;
		   ind_loop[1]=  2;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];

  
		   /// Interpolation
		   interpolation3_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd],iptstk,pos);
		   

		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];

		   pos=7;
		   /// Interpolation
		   interpolation3_(idir,param_int[nd],param_real[nd],ind_loop,iptro_p1[nd],iptro[nd],iptstk,pos); 

	
		}	     
	    }    

      

	  if ((nstep+1)%cycle==0 and param_int[nd][LEVEL] < cycle) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]<param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// Recuperation des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  param_int[nd][IJKV]-4;
		   ind_loop[1]=  param_int[nd][IJKV];
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);
	   		   
   
		}

	      else if(param_int[nd][LEVELG]>param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  //// Recuperation des flux
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=2;
		   ind_loop[0]=  1;
		   ind_loop[1]=  5;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copyflux_( idir, param_int[nd], ind_loop, iptdrodm  + shift_zone[nd], iptdrodmstk,ind,nd);
	   
		   
		}
	    }
	  
  	  
	  if ((nstep+1)%cycle==0 and param_int[nd][LEVEL]==cycle and (nstep/cycle)%2==1) /// savoir si on est a la bonne it. et si on prend les bonnes zones.
	    {

	      // cout << "coucou1" << endl;
	     
	      if(param_int[nd][LEVELG]<param_int[nd][LEVEL] and param_int[nd][LEVEL]>param_int[nd][LEVELD]) /// savoir quelles colonnes on prend
		{
		  
		  //// Stockage des valeurs
		  E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;E_Int pos=8;
		   ind_loop[0]=  1;
		   ind_loop[1]=  2;
		   ind_loop[2]=  1;
		   ind_loop[3]=  param_int[nd][IJKV+1];
		   ind_loop[4]=  1;
		   ind_loop[5]=  param_int[nd][IJKV+2];
       
		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);

		   ind=2; pos=6;

		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);
	   		   
   
		  //// Stockage des valeurs
		   pos=9;ind=1;
		   ind_loop[0]=  param_int[nd][IJKV]-1;
		   ind_loop[1]=  param_int[nd][IJKV];

		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);

		   ind=2; pos=7;

		   copy_( idir, param_int[nd], ind_loop, iptro_p1[nd], iptstk,ind,pos);
       
   
		   
		}
	    }
	  
   
  
	}
      }











 
 //////////////////////////////////// Conservativite du schema à pas de temps local d ordre 3 ////////////////////////////////////////////
  


     /* 
 

  if(param_int[0][RK]==3 and param_int[0][EXPLOC]==2 and k==2)
     {
       for (E_Int nd = 0; nd < nidom; nd++)
	 {
	   cycle=param_int[0][NSSITER]/param_int[nd][LEVEL];
	  
	   if (cycle == 4) //////////////////////////////////// Zones de plus petit niveau en temps ////////////////////////////////////////////////// 
	     {

	       //cout << nd << endl;
	       if (nstep%cycle==1 or nstep%cycle==cycle/2 or nstep%cycle==cycle-1 )
		 {

		       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;E_Int pos=0; 
		       ind_loop[2]=  1;
		       ind_loop[3]=  param_int[nd][IJKV+1];
		       ind_loop[4]=  1;
		       ind_loop[5]=  param_int[nd][IJKV+2];

		       if ((nstep/cycle)%2 == 0)
			 {
			   ind=2;
			 }
		       else
			 {
			   ind=1;
			 }
		       
		       if (param_int[nd][LEVELG]==0)  /// Stockage des flux à droite pour la conservativite
			 { 
			   ind_loop[0]= param_int[nd][IJKV]+1; ind_loop[1] =  param_int[nd][IJKV]+1;pos=0;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);
			 }

		       else if (param_int[nd][LEVELD]==0)  /// Stockage des flux à gauche pour la conservativite
			 { 

			   ind_loop[0]= 0; ind_loop[1] =  0;pos=1;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);
			 }

		       else /// Stockage des flux des 2 côtés
			 {


			   ind_loop[0]= param_int[nd][IJKV]+1; ind_loop[1] =  param_int[nd][IJKV]+1; pos=0;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);


			   //cout << ind_loop[3]  << "  " << ind_loop[5]  << endl;
			   ind_loop[0]= 0; ind_loop[1] =  0; pos=1;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);

			 }

		 }

	     }
	   else if (cycle == param_int[nd][NSSITER])  ///// Zone de plus grand pas de temps
	     {

	       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;E_Int pos=0; 
	       ind_loop[2]=  1;
	       ind_loop[3]=  param_int[nd][IJKV+1];
	       ind_loop[4]=  1;
	       ind_loop[5]=  param_int[nd][IJKV+2];

	       if (param_int[nd][LEVELG]==0) 
		 {

      
		       if (nstep%cycle==1 or nstep%cycle==cycle/2 or nstep%cycle==cycle-1)  /// Stockage des flux à droite pour la conservativite
			 { 
			   ind_loop[0]= param_int[nd][IJKV]+1; ind_loop[1] =  param_int[nd][IJKV]+1;ind=2;pos=0;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);
			 }

		       if (nstep%cycle == 0) /// Recuperation des flux de la zone de droite pour assurer la conservativite
			 {
			   //cout << "nd+1= " << nd+1 << endl;
			   //cout << ind_loop[3]  << "  " << ind_loop[5]  << endl;
			   ind_loop[0]= param_int[nd][IJKV]; ind_loop[1]=param_int[nd][IJKV];ind=2;pos=0;
			   //cout << "coe= "<< *iptcoe+shift_coe[nd] << endl; 
			   conservrk3local2_(idir,param_int[nd],param_int[nd+1],iptcoe+shift_coe[nd],iptcoe+shift_coe[nd+1],param_real[nd],ind_loop,iptro[nd],iptcstk+nd*taille,iptcstk+(nd+1)*taille,taille,nstep);

			 }
		   }


		if (param_int[nd][LEVELD]==0)  
		  {

		       if (nstep%cycle==1 or nstep%cycle==cycle/2 or nstep%cycle==cycle-1)  /// Stockage des flux à gauche pour la conservativite
			 { 
			   ind_loop[0]= 0; ind_loop[1] =  0;pos=1;ind=2;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);
			 }

		       if (nstep%cycle == 0) /// Recuperation des flux de la zone de gauche pour assurer la conservativite
			 {
			   ind_loop[0]= 1; ind_loop[1] =  1;pos=1;ind=2;			  
			   conservrk3local2_(idir,param_int[nd],param_int[nd-1],iptcoe+shift_coe[nd],iptcoe+shift_coe[nd-1],param_real[nd],ind_loop,iptro[nd],iptcstk+nd*taille,iptcstk+(nd-1)*taille,taille,nstep);

			 }


		  }



	      }

	   else   ///// Zones de  pas de temps intermédiaires
	     {

	       E_Int idir = 2; E_Int ind_loop[6]; E_Int ind=1;E_Int pos=0; 
	       ind_loop[0]=  param_int[nd][IJKV]+1;
	       ind_loop[1]=  param_int[nd][IJKV]+1;
	       ind_loop[2]=  1;
	       ind_loop[3]=  param_int[nd][IJKV+1];
	       ind_loop[4]=  1;
	       ind_loop[5]=  param_int[nd][IJKV+2];

	      if (param_int[nd][LEVELG] < param_int[nd][LEVEL] and param_int[nd][LEVEL] < param_int[nd][LEVELD])
		 {

		   if (nstep%cycle==1 or nstep%cycle==cycle/2 or nstep%cycle==cycle-1 )
			 {

			   if ((nstep/cycle)%2 == 0)
			     {
			       ind=2;
			     }
			   else
			     {
			       ind=1;
			     }
			   //cout <<"ind= " << ind << endl;
			   ///// Stockage des flux à gauche pour la conservativite
			   ind_loop[0]= 0; ind_loop[1] =  0;pos=1;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);

			   //// Stockage des flux à droite pour la conservativite
			   ind_loop[0]= param_int[nd][IJKV]+1; ind_loop[1] =  param_int[nd][IJKV]+1;pos=0;ind=2;
			   conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);			   
	  			   

			 }

		       if (nstep%cycle==0) /// Recuperation des flux de la zone de droite pour assurer la conservativite
			 {
			   ind_loop[0]= param_int[nd][IJKV]; ind_loop[1]=param_int[nd][IJKV];ind=1;pos=0;
			   conservrk3local2_(idir,param_int[nd],param_int[nd+1],iptcoe+shift_coe[nd],iptcoe+shift_coe[nd+1],param_real[nd],ind_loop,iptro[nd],iptcstk+nd*taille,iptcstk+(nd+1)*taille,taille,nstep);
			  }


		 }


	      if (param_int[nd][LEVELG] > param_int[nd][LEVEL] and param_int[nd][LEVEL] > param_int[nd][LEVELD])
		 {
		   if (nstep%cycle==1 or nstep%cycle==cycle/2 or nstep%cycle==cycle-1 )
		     {
		       if ((nstep/cycle)%2 == 0)
			 {
			   ind=2;
			 }
		       else
			 {
			   ind=1;
			 }
		       //// Stockage des flux à droite pour la conservativite
		       ind_loop[0]= param_int[nd][IJKV]+1; ind_loop[1] =  param_int[nd][IJKV]+1;pos=0;
		       conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);
		 
		       ///// Stockage des flux à gauche pour la conservativite
		       ind_loop[0]= 0; ind_loop[1] =  0;pos=1;ind=2;
		       conservrk3local_(idir,param_int[nd],ind_loop,iptdrodm + shift_zone[nd],iptcoe+shift_coe[nd], iptcstk+nd*taille,pos,taille,nstep,ind);
		     }

		   if (nstep%cycle==0) /// Recuperation des flux de la zone de gauche pour assurer la conservativite
		     {

		     ind_loop[0]= 1; ind_loop[1] =  1;pos=1;ind=1;			  
		     conservrk3local2_(idir,param_int[nd],param_int[nd-1],iptcoe+shift_coe[nd],iptcoe+shift_coe[nd-1],param_real[nd],ind_loop,iptro[nd],iptcstk+nd*taille,iptcstk+(nd-1)*taille,taille,nstep);

		     }

		   }




		 }

     

	 }
     }

 
     */













  

  //delete [] iptdrodm;
 delete [] param_int;
   

  RELEASESHAREDN( coeArray,coe);
  RELEASESHAREDN( constk  , cstk );
  RELEASESHAREDN( stock  , stk );
  RELEASESHAREDN( drodmArray,drodm);
  RELEASESHAREDN( drodmstock  , drodmstk );
  RELEASEHOOK(hook)
 

  Py_INCREF(Py_None);
  return Py_None;
  
} 
