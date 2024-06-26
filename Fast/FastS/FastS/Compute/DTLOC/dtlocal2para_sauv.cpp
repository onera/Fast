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
# include <omp.h>
using namespace std;
using namespace K_FLD;


//=============================================================================
// Idem: in place + from zone + tc compact au niveau base 
//=============================================================================
PyObject* K_FASTS::dtlocal2para_(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *work;
  PyObject *pyParam_int, *pyParam_real;
  PyObject* drodmstock;
  PyObject* constk;
  PyObject* stock;
  E_Int loc, nstep, vartype;
  E_Int omp_mode,taille_tabs, NoTransfert;

  if (!PYPARSETUPLE(args,
                    "OOOOOOOOlllll", "OOOOOOOOiiiii",
                    "OOOOOOOOlllll", "OOOOOOOOiiiii",
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work,&stock, &drodmstock, &constk,
                    &vartype, &nstep, &omp_mode, &taille_tabs, &NoTransfert))
  {
      return NULL;
  }


  vector<PyArrayObject*> hook;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);
  

    
  //E_Int NoTransfert = 3;

  //// Recuperation du tableau param_int de l'arbre t et des veceurs rop et roptmp
  E_Int**   param_intt = new E_Int*[nidomR];
  E_Float** param_realt = new E_Float*[nidomR];
  E_Float** iptro_p1 = new E_Float*[nidomR];
  E_Float** iptro    = new E_Float*[nidomR];

  for (E_Int nd = 0; nd < nidomR; nd++)

     {   
       PyObject* zone = PyList_GetItem(zonesR, nd);
 
       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
       param_intt[nd]      = K_PYTREE::getValueAI(o, hook);

	               o   = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
       param_realt[nd]     = K_PYTREE::getValueAF(o, hook);

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

  /// Recuperation du tableau de stockage des flux
  FldArrayF* drodmstk;
  K_NUMPY::getFromNumpyArray(drodmstock,drodmstk, true); E_Float* iptdrodmstk = drodmstk->begin();

  /// Recuperation du tableau de stockage des flux pour conservativite
  FldArrayF* cstk;
  K_NUMPY::getFromNumpyArray(constk, cstk, true); E_Float* iptcstk = cstk->begin();

  /// Tableau de travail communs explicite/implicite
  PyObject* drodmArray = PyDict_GetItemString(work,"rhs"); FldArrayF* drodm;
  K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();
  
  // Tableau de travail coe   ( dt/vol et diags LU)
  PyObject* coeArray = PyDict_GetItemString(work,"coe"); FldArrayF* coe;
  K_NUMPY::getFromNumpyArray(coeArray, coe, true); E_Float* iptcoe = coe->begin();

  /*-------------------------------------*/
  /* Extraction tableau int et real      */
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();

  E_Int nvars;

  if( vartype <= 3 &&  vartype >= 1) nvars =5;
  else                               nvars =6;


  //cout << "ipt_param_int[0]= " << ipt_param_int[0] << endl;
  ///// Les shifts pour les zones /////

  E_Int nbcomIBC = ipt_param_int[1];
  E_Int nbcomID  = ipt_param_int[2+nbcomIBC];
  
  E_Int shift_graph = nbcomIBC + nbcomID + 2;

  //E_Int threadmax_sdm  = __NUMTHREADS__;
  E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];

 
  E_Int nrac = ipt_param_int[ ech +1 ];


  E_Int a=0;
  E_Int b=0;
  E_Int shift_zone[nidomR];
  E_Int shift_coe [nidomR];

     for (E_Int nd = 0; nd < nidomR; nd++)
       {
	 shift_zone[nd]=a;
	 a=a+param_intt[nd][ NDIMDX ]*param_intt[nd][ NEQ ];	 
       }
      for (E_Int nd = 0; nd < nidomR; nd++)
       {
	 shift_coe[nd]=b;
	 b=b+param_intt[nd][ NDIMDX ]*param_intt[nd][ NEQ_COE ];	 
       }

      E_Int taille=taille_tabs/nrac;


  vector<E_Int> nbraczone(nidomR,0);
  vector<E_Int> stockzone(nidomR,-1);


#pragma omp parallel default(shared) //private(cycle)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
   E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;
   if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
   else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

   E_Int topology[3];
   topology[0]=0;
   topology[1]=0;
   topology[2]=0;

  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
    {

      E_Int timelevel = ipt_param_int[ ech +3 ]; 
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];
      
      E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 26*irac;

      E_Int cycle;
      //cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];
      cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];

      //cout << irac <<"  "<< ipt_param_int[debut_rac + 25] <<"  "<< ipt_param_int[debut_rac + 24]  << endl;

      //if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
      if (ipt_param_int[debut_rac + 25] > ipt_param_int[debut_rac + 24])

       {

	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 

	     E_Int donorPts_[6]; 
	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
             donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
             donorPts_[2] =  ipt_param_int[debut_rac + 9];
             donorPts_[3] =  ipt_param_int[debut_rac + 10];
             donorPts_[4] =  ipt_param_int[debut_rac + 11];
             donorPts_[5] =  ipt_param_int[debut_rac + 12];
	     int dir = ipt_param_int[debut_rac + 13];
	     E_Int profondeur = ipt_param_int[debut_rac + 20];

	     //cout << dir << endl;

	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);

	     E_Int taillefenetre;
	     taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);

	     E_Int ipt_ind_dm_omp_thread[6]; 

	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                              donorPts_,
                              topology, ipt_ind_dm_omp_thread);

	     E_Int donorPts[6];
	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];

	     //cout << "pas de temps zone donneuse + petit  "<< NoD <<" "<< NoTransfert << endl;

	     if (nstep%cycle==1 and (nstep/cycle)%2==0) /// La 1ere sous-iteration du cycle
	       {

		 //cout <<donorPts[0] <<" "<< donorPts[1]<< " "<<  donorPts[2]<<" "<< donorPts[3]<<" "<< donorPts[4]<<" "<< donorPts[5]<<" "<<NoD<< endl;
		 
		 //cout << ech + 4 + timelevel*2 + 1 + nrac*16 + 21*irac << endl;
		 
		 //// stockage des flux en vue de l'interpolation pour la zone de + gd pas de temps (ici NoR) en position 0 dans le tableau de stockage des flux (raccord)
		 E_Int ind=1;
		 copyflux_rk3localpara_(param_intt[NoD],donorPts, donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);

		 /// stockage de yn en position 0 dans le tableau de stockage du raccord
		 copy_rk3localpara_(param_intt[NoD], donorPts,donorPts_, iptro[NoD], iptstk + irac*taille + 0*taillefenetre,ind,taillefenetre); 



	       }
	 
	     if (nstep%cycle==cycle/2 and (nstep/cycle)%2==0 ) /// Recuperation des valeurs stockées en position 0 (yn ou yinterpolé suivant la parité du cycle) et stockage de y2 ou y6
		 
	       {   

		 //// Recuperation et transformation de f(yn) + coeff*f(y1) pour obtenir alpha*f(yn) + beta*f(y1)
		 E_Int ind=2;
		 copyflux_rk3localpara_(param_intt[NoD], donorPts , donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);
			 
	       }


	     if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)

	       {

		 /// Stockage de y6  en position 1 dans le tableau de stockage du raccord
		 E_Int ind=1;
		 copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, iptro[NoD], iptstk + irac*taille + 1*5*taillefenetre,ind,taillefenetre); 
 
	       }


	     // if (nstep%cycle==0) /// Interpolation pour la zone adjacente de pas de temps plus petit et switch pointeurs
	     // {  
		   
	     // }


       }

      //else if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
      else if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24])
       {

	  E_Int donorPts_[6];  E_Int idir = 2;
	  donorPts_[0] =  ipt_param_int[debut_rac + 7];
          donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
          donorPts_[2] =  ipt_param_int[debut_rac + 9];
          donorPts_[3] =  ipt_param_int[debut_rac + 10];
          donorPts_[4] =  ipt_param_int[debut_rac + 11];
          donorPts_[5] =  ipt_param_int[debut_rac + 12];
	  int dir = ipt_param_int[debut_rac + 13];
	  E_Int profondeur = ipt_param_int[debut_rac + 20];

	  donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);
		       
	  E_Int taillefenetre;
	  taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);

	  E_Int ipt_ind_dm_omp_thread[6]; 

	  indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                              donorPts_,
                              topology, ipt_ind_dm_omp_thread);

	  E_Int donorPts[6];
	  donorPts[0]=ipt_ind_dm_omp_thread[0];
	  donorPts[1]=ipt_ind_dm_omp_thread[1];
	  donorPts[2]=ipt_ind_dm_omp_thread[2];
	  donorPts[3]=ipt_ind_dm_omp_thread[3];
	  donorPts[4]=ipt_ind_dm_omp_thread[4];
	  donorPts[5]=ipt_ind_dm_omp_thread[5];

	  //cout << "pas de temps zone donneuse + grand  "<< NoD <<" "<< NoTransfert << endl;
	     
	  if (nstep%cycle == 1)

	    {

	      //cout << donorPts[0]<<" "<< donorPts[1]<<"  "<< donorPts[2]<<"  "<< donorPts[3]<<"  "<< donorPts[4]<<"  "<< donorPts[5]<<" "<<"  "<< NoD << endl; 

	      //cout << ech + 4 + timelevel*2 + 1 + nrac*16 + 21*irac << endl;

	      /// Stockage de yn en position 0 dans le tableau de stocakge du raccord 
	      E_Int ind=1;
	      copy_rk3localpara_(param_intt[NoD], donorPts , donorPts_, iptro[NoD], iptstk + irac*taille + 0*taillefenetre,ind,taillefenetre); 

	      /// Stockage de y3 en position 1 dans le tableau de stockage du raccord 
	      copy_rk3localpara_(param_intt[NoD], donorPts , donorPts_, iptro_p1[NoD], iptstk + irac*taille + 1*5*taillefenetre,ind,taillefenetre); 
    
	      /// stockage de drodm (yn)
	      ind=1;
	      copyflux_rk3local2para_(param_intt[NoD], donorPts , donorPts_ ,iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);

	    }
	      
	      
	  if (nstep%cycle == cycle/4)

	    {		       
	      
	      /// Interpolation de y2
	      interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts, donorPts_, iptstk + irac*taille + 0*taillefenetre,iptstk + irac*taille + 1*5*taillefenetre, iptro[NoD],dir,taillefenetre,nstep);
	    }

	      
	     
	  if (nstep%cycle == cycle/2-1)

	    {
	       
	      /// Recuperation de y3 en position 1 dans le tableau de stockage du raccord 
	      E_Int ind=2;
	      copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, iptro_p1[NoD], iptstk + irac*taille + 1*5*taillefenetre,ind,taillefenetre); 

	    }


	  if (nstep%cycle == cycle/2)

	    {

	      /// Switch pointeurs

	      /// Stockage de y4 en position 4 dans le tableau de stockage du raccord 
	      E_Int ind=1;
	      copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, iptro[NoD], iptstk + irac*taille + 2*5*taillefenetre,ind,taillefenetre);

	      /// stockage de drodm(y1) pour obtenir alpha*drodm(yn) + beta*drodm(y1)
	      ind=2;
	      copyflux_rk3local2para_(param_intt[NoD], donorPts , donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);


	    }

	  if (nstep%cycle == cycle/2 + 1)

	    {
	      /// Switch pointeurs

	      //Interpolation de y6
	      interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre,iptstk + irac*taille + 1*5*taillefenetre, iptro_p1[NoD],dir,taillefenetre,nstep);

	    }


	  // if (nstep%cycle == 0) 
	  // {

	  /// Switch pointeurs

	  // }
	    

       }
     //else

     //{

     // if (nstep%cycle==0) /// Interpolation pour la zone adjacente de pas de temps plus petit et switch pointeurs
     //  {


	     // if (stockzone[NoD] != NoD and nstep != param_intt[NoD][NSSITER])
	     //	 {
		 

	     //   stockzone[NoD] = NoD;
	     //	 }


     //  }
	        
     // }
	  

    }// boucle raccords

 }// fin zone omp 

	 


#pragma omp parallel default(shared) //private(cycle)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
   E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;
   if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
   else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

   E_Int topology[3];
   topology[0]=0;
   topology[1]=0;
   topology[2]=0;



  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

    {

      //E_Int shift_rac =  ech + 2 + irac;
      E_Int timelevel = ipt_param_int[ ech +3 ]; 
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];


      vector<E_Float*> vectOfRcvFields(nvars);
      vector<E_Float*> vectOfDnrFields(nvars);

      for (E_Int eq = 0; eq < nvars; eq++)
	{
	  vectOfRcvFields[eq] = iptro[ NoR] + eq*param_intt[ NoR ][NDIMDX];
	  vectOfDnrFields[eq] = iptro[ NoD] + eq*param_intt[ NoD ][NDIMDX];
	}


      //E_Int cycle;  
      //cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];

      E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 26*irac;

      E_Int cycle;
      cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];

      //if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
      if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24])
       {

	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 

	     E_Int donorPts_[6]; 
	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
             donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
             donorPts_[2] =  ipt_param_int[debut_rac + 9];
             donorPts_[3] =  ipt_param_int[debut_rac + 10];
             donorPts_[4] =  ipt_param_int[debut_rac + 11];
             donorPts_[5] =  ipt_param_int[debut_rac + 12];
	     int dir = ipt_param_int[debut_rac + 13];
	     E_Int profondeur = ipt_param_int[debut_rac + 20];

	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);

	     //cout << "irac= " <<irac<< endl;
	     E_Int taillefenetre;
	     taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);

	     E_Int ipt_ind_dm_omp_thread[6]; 

	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                             donorPts_,
                             topology, ipt_ind_dm_omp_thread);

	     E_Int donorPts[6];
	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];


	     if (nstep%cycle == 1)

	       {

		 /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn 
		  interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre,iptstk + irac*taille + 1*5*taillefenetre, iptro_p1[NoD],dir,taillefenetre,nstep);
	     
	       }
	     

	     if (nstep%cycle == cycle/2 + cycle/4)
 
	       {

		 //Interpolation de y6 ds rop
		 E_Float coeff=1.0;
		 interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre, iptdrodmstk + irac*taille + 0*taillefenetre, iptro[NoD],taillefenetre,coeff);


		       
#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

		 /* 
       E_Int pos;
       pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
       pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
       pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts   = ipt_param_int +  pos;
       pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
       pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

       E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ]; E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;

       E_Int ideb      = 0;
       E_Int ifin      = 0;
       E_Int shiftCoef = 0;
       E_Int shiftDonnor  = 0;
       E_Int sizecoefs = 0;
       
       for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++)
       { 
        E_Int type      = types[ifin];

        //SIZECF(type, meshtype, sizecoefs);
        ifin =  ifin + ntype[ 1 + ndtyp];

        E_Int pt_deb, pt_fin;

        // Calcul du nombre de champs a traiter par chaque thread
        E_Int size_bc =  ifin-ideb;
        E_Int chunk   =  size_bc/Nbre_thread_actif;
        E_Int r       =  size_bc - chunk*Nbre_thread_actif;
        // pts traitees par thread
        if (ithread <= r)
             { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
        else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; } 

        //Si type 0, calcul sequentiel
        if      ( type == 0 )
          { if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
            else             { pt_deb = ideb; pt_fin = ideb;}
          }

          E_Int noi       = shiftDonnor;                             // compteur sur le tableau d indices donneur
          E_Int indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;

#         include "commonTransfert_5eq.h" 

  //        } //chunk
  //
          ideb       = ideb + ifin;
          shiftCoef  = shiftCoef    +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
          shiftDonnor= shiftDonnor  +  ntype[1+ndtyp];           //shift donnor entre 2 types successif
       }// type 

		 */



		 /*

		 E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];
		 E_Int pos;
		 //pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
		 //pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
		 pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts2   = ipt_param_int +  pos;
		 pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
		 //pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

		 E_Int size_bc =  nbRcvPts ;
		 E_Int chunk   =  size_bc/Nbre_thread_actif;
		 E_Int r       =  size_bc - chunk*Nbre_thread_actif;
		 E_Int ideb    =  0;
		 E_Int pt_deb;
		 E_Int pt_fin;
		 // pts traitees par thread
		 if (Nbre_thread_actif != 1)
		   {
		     if (ithread <= r)
		       { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); pt_fin=pt_fin-1;}
		     else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk;  pt_fin=pt_fin-1;} 

		   }
		 else
		   {
		     pt_deb=0;
		     pt_fin=nbRcvPts-1;
		     nbRcvPts = pt_fin - pt_deb +1;

		   }

		 remp_cellfictivespara_(param_intt[NoD],param_intt[NoR], iptro[NoD], iptro[NoR],donorPts2,rcvPts,nbRcvPts,pt_deb,pt_fin);
		 */





	       }


      }

      //else if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 
     else if (ipt_param_int[debut_rac + 25] > ipt_param_int[debut_rac + 24])
       {

	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	     E_Int donorPts_[6]; 
	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
             donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
             donorPts_[2] =  ipt_param_int[debut_rac + 9];
             donorPts_[3] =  ipt_param_int[debut_rac + 10];
             donorPts_[4] =  ipt_param_int[debut_rac + 11];
             donorPts_[5] =  ipt_param_int[debut_rac + 12];
	     int dir = ipt_param_int[debut_rac + 13];
	     E_Int profondeur = ipt_param_int[debut_rac + 20];
	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);

	     //cout << "irac= " <<irac<< endl;
	      E_Int taillefenetre;
	      taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);

	     E_Int ipt_ind_dm_omp_thread[6]; 

	      indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                               donorPts_,
                               topology, ipt_ind_dm_omp_thread);
	      E_Int donorPts[6];
	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];


	     if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)

	       {


		 /// Ecriture de yinterpolé dans rop
		 E_Float coeff=2.0 ;		        	       
		 interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre, iptdrodmstk + irac*taille + 0*taillefenetre,iptro[NoD],taillefenetre,coeff);  


#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]
		 /*

       E_Int pos;
       pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
       pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
       pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts   = ipt_param_int +  pos;
       pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
       pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

       E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ]; E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;
       //if(ibc ==1)
       // { 
       //  xPC     = ptrCoefs + nbInterpD;
       // xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
       // xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
       //  densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
       // }

       E_Int ideb      = 0;
       E_Int ifin      = 0;
       E_Int shiftCoef = 0;
       E_Int shiftDonnor  = 0;
       E_Int sizecoefs = 0;
       
       for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++)
       { 
        E_Int type      = types[ifin];

        //SIZECF(type, meshtype, sizecoefs);
        ifin =  ifin + ntype[ 1 + ndtyp];

        E_Int pt_deb, pt_fin;

        // Calcul du nombre de champs a traiter par chaque thread
        E_Int size_bc =  ifin-ideb;
        E_Int chunk   =  size_bc/Nbre_thread_actif;
        E_Int r       =  size_bc - chunk*Nbre_thread_actif;
        // pts traitees par thread
        if (ithread <= r)
             { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
        else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; } 

        //Si type 0, calcul sequentiel
        if      ( type == 0 )
          { if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
            else             { pt_deb = ideb; pt_fin = ideb;}
          }

          E_Int noi       = shiftDonnor;                             // compteur sur le tableau d indices donneur
          E_Int indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;

#         include "commonTransfert_5eq.h" 

  //        } //chunk
  //
          ideb       = ideb + ifin;
          shiftCoef  = shiftCoef    +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
          shiftDonnor= shiftDonnor  +  ntype[1+ndtyp];           //shift donnor entre 2 types successif
       }// type 

		 */

		 /*
		 E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];
	     
		 pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
		 pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
		 pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts2   = ipt_param_int +  pos;
		 pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
		 // //pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

		 // //cout << irac << endl;

		 E_Int size_bc =  nbRcvPts ;
		 E_Int chunk   =  size_bc/Nbre_thread_actif;
		 E_Int r       =  size_bc - chunk*Nbre_thread_actif;
		 E_Int ideb    =  0;
		 E_Int pt_deb;
		 E_Int pt_fin;
		 //pts traitees par thread
		 if (Nbre_thread_actif != 1)
		   {
		     if (ithread <= r)
		       { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); pt_fin=pt_fin-1;}
		     else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk;  pt_fin=pt_fin-1;} 
		   }
		 else
		   {
		     pt_deb=0;
		     pt_fin=nbRcvPts-1;
		     nbRcvPts = pt_fin - pt_deb +1;
		   }

	     
		 remp_cellfictivespara_(param_intt[NoD],param_intt[NoR], iptro[NoD], iptro[NoR],donorPts2,rcvPts,nbRcvPts,pt_deb,pt_fin);

		 */
	       }
	     
       }

	      

    }

  } // fin zone omp



  //E_Int taille_=2000000/nrac;

  

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////// Ajout de la conservativite ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
 

#pragma omp parallel default(shared) //private(cycle)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
   E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int Nbre_thread_actif_loc, ithread_loc;
   if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
   else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

   E_Int topology[3];
   topology[0]=0;
   topology[1]=0;
   topology[2]=0;


  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

    {

      E_Int timelevel = ipt_param_int[ ech +3 ]; 
      //E_Int shift_rac =  ech + 2 + irac;
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];

      E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 26*irac;       

      //E_Int cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];
      //E_Int cycleR = param_intt[NoR][NSSITER]/param_intt[NoR][LEVEL];

      E_Int cycle  = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];
      E_Int cycleR = param_intt[NoR][NSSITER]/ipt_param_int[debut_rac + 24];


      //  if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
      if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24]) 
       {

	     E_Int pos;E_Int ind;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	     E_Int donorPts_[6]; 
	     donorPts_[0] =  ipt_param_int[debut_rac + 0];
             donorPts_[1] =  ipt_param_int[debut_rac + 1];
             donorPts_[2] =  ipt_param_int[debut_rac + 2];
             donorPts_[3] =  ipt_param_int[debut_rac + 3];
             donorPts_[4] =  ipt_param_int[debut_rac + 4];
             donorPts_[5] =  ipt_param_int[debut_rac + 5];
	     int dir = ipt_param_int[debut_rac + 6];


	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

	     //cout << "receveurs_cons= " << donorPts_[0]<<" "<< donorPts_[1] <<" "<< donorPts_[2]<<" "<< donorPts_[3] << endl;

	     E_Int taillefenetreR;
	     taillefenetreR = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 

	     E_Int ipt_ind_dm_omp_thread[6]; 

	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                              donorPts_,
                              topology, ipt_ind_dm_omp_thread);
	     E_Int donorPts[6];
	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];


	     if (nstep%cycleR==1 or nstep%cycleR==cycleR/2 or nstep%cycleR==cycleR-1)
	       {


		 if ((nstep/cycleR)%2 == 0)
		   {
		      ind=2;
		   }
		 else
		   {
		      ind=1;
		   }
		 //cout << "receveurs= " << donorPts_[0]<<" "<< donorPts_[1] <<" "<< donorPts_[2]<<" "<< donorPts_[3] << endl;
		 //cout << "ind= " << ind << endl;
		 conservrk3local3para_(param_intt[NoR],donorPts,donorPts_,iptdrodm + shift_zone[NoR],iptcoe+shift_coe[NoR], iptcstk+irac*taille,taillefenetreR,nstep,ind);

	       }


	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	      


	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
             donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
             donorPts_[2] =  ipt_param_int[debut_rac + 9];
             donorPts_[3] =  ipt_param_int[debut_rac + 10];
             donorPts_[4] =  ipt_param_int[debut_rac + 11];
             donorPts_[5] =  ipt_param_int[debut_rac + 12];
	     dir = ipt_param_int[debut_rac + 13];

	     //cout << "donneurs= "<<  donorPts_[0] <<" "<< donorPts_[1] <<" "<< donorPts_[2] <<" "<< donorPts_[3] <<" "<<dir<< endl;


	     //donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_intt[NoD][NIJK]-1-4) ;
	     //donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_intt[NoD][NIJK]-1-4);


	     
	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

	     //cout << "donneurs_conserv= "<<  donorPts_[0] <<" "<< donorPts_[1] <<" "<< donorPts_[2] <<" "<< donorPts_[3] << endl;

	     E_Int taillefenetreD;	     
	     taillefenetreD = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 

	     
	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                              donorPts_,
                              topology, ipt_ind_dm_omp_thread);

	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];

	     if (nstep%cycle==1 or nstep%cycle==cycle/2 or nstep%cycle==cycle-1)
	       {


		E_Int ind=2;
		conservrk3local3para_(param_intt[NoD],donorPts,donorPts_,iptdrodm + shift_zone[NoD],iptcoe+shift_coe[NoD], iptcstk+irac*taille+5*taillefenetreR,taillefenetreD,nstep,ind);

	       }


	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
             donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
             donorPts_[2] =  ipt_param_int[debut_rac + 9];
             donorPts_[3] =  ipt_param_int[debut_rac + 10];
             donorPts_[4] =  ipt_param_int[debut_rac + 11];
             donorPts_[5] =  ipt_param_int[debut_rac + 12];
	     //int dir = ipt_param_int[ech + 4 + timelevel*2 + 1 + nrac*16 + 14*irac + 13];

	     //donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_intt[NoD][NIJK]-4-1) ;
	     //donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_intt[NoD][NIJK]-4-1);

	     taillefenetreD = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 


	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                              donorPts_,
                              topology, ipt_ind_dm_omp_thread);

	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];

	     //dir = ipt_param_int[ech + 4 + timelevel*2 + 1 + nrac*16 + 20*irac + 13];
	     //donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     //donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];


	     E_Int donorPts__[6]; 
	     donorPts__[0] =  ipt_param_int[debut_rac + 0];
             donorPts__[1] =  ipt_param_int[debut_rac + 1] ;
             donorPts__[2] =  ipt_param_int[debut_rac + 2];
             donorPts__[3] =  ipt_param_int[debut_rac + 3];
             donorPts__[4] =  ipt_param_int[debut_rac + 4];
             donorPts__[5] =  ipt_param_int[debut_rac + 5];
	     //dir = ipt_param_int[ech + 4 +timelevel*2 + 1 + nrac*16 + 20*irac + 6];
	     //donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     //donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

	     E_Int* transfo;
	     transfo = &ipt_param_int[debut_rac + 14];
  

             E_Int trans1 =  ipt_param_int[debut_rac + 14];
             E_Int trans2 =  ipt_param_int[debut_rac + 15];
	     E_Int trans3 =  ipt_param_int[debut_rac + 16];
	       

             E_Int pt1 =  ipt_param_int[debut_rac + 17];
	     E_Int pt2 =  ipt_param_int[debut_rac + 18];
	     E_Int pt3 =  ipt_param_int[debut_rac + 19];

             E_Int ratio1 =  ipt_param_int[debut_rac + 21];
	     E_Int ratio2 =  ipt_param_int[debut_rac + 22];
	     E_Int ratio3 =  ipt_param_int[debut_rac + 23];

	     //cout << ratio1 <<" "<<ratio2<<" "<<ratio3<< endl;


	     if (nstep%cycle==cycle-1) /// Opération pour assurer la conservativité
	       {


#pragma omp barrier // On attend que tous les threads aient calculé leur bilan de flux
		 //cout << "zone "<< NoD<< endl;

		 //cout <<  "receveurs= " <<  donorPts__[0]<<" "<< donorPts__[1] <<" "<< donorPts__[2]<<" "<< donorPts__[3] <<" "<< donorPts__[4]<<" "<< donorPts__[5]<<endl;

		 //cout << "transfo= "    << *transfo << " "<<*(transfo+1)<<" "<<*(transfo+2) << " "<< irac << endl;
		 //cout << NoD  <<" "<<  NoR  << " "<< irac << endl;
		 //cout << "dir= " << dir << endl;
		 //cout << "donneurs= "   <<  donorPts_[0]<<" "<< donorPts_[1] <<" "<< donorPts_[2]<<" "<< donorPts_[3] <<" "<< donorPts_[4]<<" "<< donorPts_[5]<<endl;

		 conservrk3local4para_(param_intt[NoD],param_intt[NoR],iptcoe+shift_coe[NoD],iptcoe+shift_coe[NoR],param_realt[NoD],donorPts,donorPts_,donorPts__,iptro_p1[NoD],iptcstk+irac*taille+5*taillefenetreR,iptcstk+irac*taille,taillefenetreR,taillefenetreD,nstep,dir,pt1,pt2,pt3,transfo,ratio1,ratio2,ratio3);

	       }

       }

  


    }


}//fin zone omp
 

 RELEASESHAREDN( constk  , cstk );
 RELEASESHAREDN( stock  , stk );
 RELEASESHAREDN( drodmstock  , drodmstk );
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );
 RELEASESHAREDN( coeArray,coe);
 RELEASESHAREDN( drodmArray,drodm);

 Py_INCREF(Py_None);
 return Py_None;
}
