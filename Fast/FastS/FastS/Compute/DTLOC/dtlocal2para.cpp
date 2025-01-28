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
PyObject* K_FASTS::dtlocal2para(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *work;
  PyObject *pyParam_int, *pyParam_real; // param de tc
  PyObject* drodmstock;
  PyObject* constk;
  PyObject* stock;
  E_Int loc, nstep, vartype;
  E_Int omp_mode,taille_tabs;

  if (!PYPARSETUPLE(args,
                    "OOOOOOOOllll", "OOOOOOOOiiii",
                    "OOOOOOOOllll", "OOOOOOOOiiii",
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work,&stock, &drodmstock, &constk,
                    &vartype,&nstep,&omp_mode,&taille_tabs))
  {
      return NULL;
  }


  vector<PyArrayObject*> hook;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);
  

    
  E_Int NoTransfert = 1;

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
  //// Recuperation du tableau drodm
  PyObject* drodmArray = PyList_GetItem(work,2); FldArrayF* drodm;
  K_NUMPY::getFromNumpyArray(drodmArray, drodm, true); E_Float* iptdrodm = drodm->begin();

  /// Tableau de travail coe   ( dt/vol et diags LU)
  PyObject* coeArray = PyList_GetItem(work,1); FldArrayF* coe;
  K_NUMPY::getFromNumpyArray(coeArray, coe, true); E_Float* iptcoe = coe->begin();


 
  /*-------------------------------------*/
  /* Extraction tableau int et real de tc*/
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


  vector<E_Float*> vectOfRcvFields_rop(nvars);
  vector<E_Float*> vectOfDnrFields_rop(nvars);

  vector<E_Float*> vectOfRcvFields_roptmp(nvars);
  vector<E_Float*> vectOfDnrFields_roptmp(nvars);


  ///// Les shifts pour les zones /////
  
  E_Int ech  = ipt_param_int[ NoTransfert ];
  E_Int nrac = ipt_param_int[ ech +1 ];

  omp_mode = 0;


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

   //vector<E_Int> stockzone(nidomR,-1);  

  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
    {

      //E_Int shift_rac =  ech + 2 + irac;
      E_Int timelevel = ipt_param_int[ ech +3 ]; 
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5  ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 ];


      E_Int cycle;
      cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];

      //  if (ithread==2)
      //     {
      //       cout << cycle <<" "<<param_intt[NoD][LEVEL]<<" "<<param_intt[NoR][LEVEL] <<"  "<< irac << endl;  
      //     }

      //cout << NoD <<"  "<< NoR << endl;

     if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 
       {



	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	     E_Int donorPts_[6]; 

             E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
	     donorPts_[0] = ipt_param_int[adr + 7];
             donorPts_[1] = ipt_param_int[adr + 8] ;
             donorPts_[2] = ipt_param_int[adr + 9];
             donorPts_[3] = ipt_param_int[adr +10];
             donorPts_[4] = ipt_param_int[adr +11];
             donorPts_[5] = ipt_param_int[adr +12];
	     E_Int dir    = ipt_param_int[adr +13];
	     //cout << dir << endl;
	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

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




     
	 if (nstep%cycle==1 and (nstep/cycle)%2==0) /// La 1ere sous-iteration du cycle
	   {


	     //// stockage des flux en vue de l'interpolation pour la zone de + gd pas de temps (ici NoR) en position 0 dans le tableau de stockage des flux (raccord)
	     E_Int ind=1;
	     copyflux_rk3localpara_(param_intt[NoD],donorPts, donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);

	     /// stockage de yn en position 0 dans le tableau de stockage du raccord
	     copy_rk3localpara_(param_intt[NoD], donorPts,donorPts_, iptro[NoD], iptstk + irac*taille + 0*taillefenetre,ind,taillefenetre); 

	     cout <<donorPts[0] <<" "<< donorPts[1]<< " "<<  donorPts[2]<<" "<< donorPts[3]<<" "<< donorPts[4]<<" "<< donorPts[5]<< endl;


	   }
	 
	 if (nstep%cycle==cycle/2 and (nstep/cycle)%2==0 ) /// Recuperation des valeurs stockées en position 0 (yn ou yinterpolé suivant la parité du cycle) et stockage de y2 ou y6
		 
	   {   

	     //// Recuperation et transformation de f(yn) + coeff*f(y1) pour obtenir alpha*f(yn) + beta*f(y1)
	     E_Int ind=2;
	     copyflux_rk3localpara_(param_intt[NoD], donorPts , donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);
	     //cout << donorPts[0]<<" "<< donorPts[1]<<"  "<< donorPts[2]<<"  "<< donorPts[3]<<"  "<< donorPts[4]<<"  "<< donorPts[5]<<" "<<NoD<< endl;		       
			 
	   }


	 if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
	 //if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1)

	   {


	   /// Stockage de y6  en position 1 dans le tableau de stockage du raccord
	   E_Int ind=1;
	   copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, iptro[NoD], iptstk + irac*taille + 1*5*taillefenetre,ind,taillefenetre); 

	   /// Recuperation de yinterpolé en position 0 dans le tableau de stockage du raccord
	   //ind=2;    
	   //copy_rk3local_(param_intt[NoD], donorPts , iptro[NoD], iptstk + irac*taille + 0*taillefenetre,ind,taillefenetre); 

	   }


	 // if (nstep%cycle==0) /// Interpolation pour la zone adjacente de pas de temps plus petit et switch pointeurs
	 // {  
		   
	 //   if ((nstep/cycle)%2==1 )
	 //	     {
	 //	       E_Float coeff=2.0 ;//cout << NoD << endl;		        	       
	 //	       interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre, iptdrodmstk + irac*taille + 0*taillefenetre,taillefenetre,coeff);  
	 //	     }
	 //	 }
	       //#pragma omp barrier
       // if (ithread==2)
       //  {
       //     cout << cycle <<" "<<param_intt[NoD][LEVEL]<<" "<<param_intt[NoR][LEVEL] << " "<< irac << endl;  
       //    }


       }
 
     else if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
       {

		 //cout << "irac= " <<irac<< endl;

	 //cout << NoD << endl;
	 
	 E_Int donorPts_[6];  E_Int idir = 2;
         E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
         donorPts_[0] = ipt_param_int[adr + 7];
         donorPts_[1] = ipt_param_int[adr + 8] ;
         donorPts_[2] = ipt_param_int[adr + 9];
         donorPts_[3] = ipt_param_int[adr +10];
         donorPts_[4] = ipt_param_int[adr +11];
         donorPts_[5] = ipt_param_int[adr +12];
         E_Int dir    = ipt_param_int[adr +13];

	 donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;
		       
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

	 //if (ithread==2)
	 //    {
	 //	 cout << donorPts[0]<<" "<<donorPts[1]<<" "<<donorPts[2]<<" "<<donorPts[3]<<" "<<donorPts[4]<<" "<< donorPts[5]<<endl;  
	 //      }

	 //cout <<donorPts[0] <<" "<< donorPts[1]<< " "<<  donorPts[2]<<" "<< donorPts[3]<<" "<< donorPts[4]<<" "<< donorPts[5]<< endl;

	     
	      if (nstep%cycle == 1)

		     {

		       //cout << NoD <<" "<< cycle<<endl;

		       cout << donorPts[0]<<" "<< donorPts[1]<<"  "<< donorPts[2]<<"  "<< donorPts[3]<<"  "<< donorPts[4]<<"  "<< donorPts[5]<<" "<<NoD<< endl; 

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

		       //cout << NoD << endl;
		       /// Switch pointeurs

		       //if (stockzone[NoD] != NoD)
		       //{


			    // PyObject* zone = PyList_GetItem(zonesR, NoD);

			    // PyObject* o  = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
			    // PyObject* tp1  = K_PYTREE::getNodeFromName1( o        , "Density_P1");

			    // o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
			    // PyObject* t1  = K_PYTREE::getNodeFromName1( o        , "Density");


			    // PyObject* ts1 = PyList_GetItem(t1,1);

			    // PyList_SetItem(t1,1, PyList_GetItem(tp1,1));


			    // PyList_SetItem(tp1,1,ts1);


			   //switchvectors_(iptro[NoD],iptro_p1[NoD],param_intt[NoD]) ;
			   //cout<< "apres switchvector" << endl;

			   //PyObject* zone = PyList_GetItem(zonesR, NoD);
 
		       // stockzone[NoD] = NoD;
		       // }


		       /// Stockage de y4 en position 4 dans le tableau de stockage du raccord 
		       E_Int ind=1;
		       //copy_rk3local_(param_intt[NoD], donorPts , iptro_p1[NoD], iptstk + irac*taille + 2*5*taillefenetre,ind,taillefenetre);

		       //cout << donorPts[0] <<" "<<  donorPts[1] <<" "<<donorPts[2]<<" "<<donorPts[3] << endl;
		       copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, iptro[NoD], iptstk + irac*taille + 2*5*taillefenetre,ind,taillefenetre);

     		       /// stockage de drodm(y1) pour obtenir alpha*drodm(yn) + beta*drodm(y1)
		       ind=2;
		       copyflux_rk3local2para_(param_intt[NoD], donorPts , donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + irac*taille + 0*taillefenetre,ind,taillefenetre);


		     }

		  if (nstep%cycle == cycle/2 + 1)
		     {

		       interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre,iptstk + irac*taille + 1*5*taillefenetre, iptro_p1[NoD],dir,taillefenetre,nstep);

		     }

		   if (nstep%cycle == cycle/2 + cycle/4) 
		     {


		       E_Float coeff=1.0;
		       interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre, iptdrodmstk + irac*taille + 0*taillefenetre, iptro[NoD],taillefenetre,coeff);


		     }

		   // if (nstep%cycle == 0) 
		   // {

		       /// Switch pointeurs
		      // if (stockzone[NoD] != NoD and nstep != param_intt[NoD][NSSITER])
		      // {
		      //switchvectors_(iptro[NoD],iptro_p1[NoD],param_intt[NoD]) ;

		      //    PyObject* zone = PyList_GetItem(zonesR, NoD);

		      //    PyObject* o  = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
		      //    PyObject* tp  = K_PYTREE::getNodeFromName1( o        , "Density_P1");

		      //    o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
		      //    PyObject* t  = K_PYTREE::getNodeFromName1( o        , "Density");

			    //cout << "pointeur density avant= " << PyList_GetItem(t,1)  << endl;

		      //    PyObject* t1 = PyList_GetItem(t,1);
		      //    PyList_SetItem(t,1, PyList_GetItem(tp,1));
		      //    PyList_SetItem(tp,1,t1);


		      // stockzone[NoD] = NoD;
		      // }
	      
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

  // }// fin zone omp 

	 


  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

    {

      //E_Int shift_rac =  ech + 2 + irac;
      E_Int timelevel = ipt_param_int[ ech +3 ]; 
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];

      E_Int cycle;  
      cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];

     if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
      {

	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	     E_Int donorPts_[6]; 
             E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
             donorPts_[0] = ipt_param_int[adr + 7];
             donorPts_[1] = ipt_param_int[adr + 8] ;
             donorPts_[2] = ipt_param_int[adr + 9];
             donorPts_[3] = ipt_param_int[adr +10];
             donorPts_[4] = ipt_param_int[adr +11];
             donorPts_[5] = ipt_param_int[adr +12];
             E_Int dir    = ipt_param_int[adr +13];
	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

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


	     //  donorPts[0]=donorPts_[0];
	     //  donorPts[1]=donorPts_[1];
	     //  donorPts[2]=donorPts_[2];
	     // donorPts[3]=donorPts_[3];
	     //  donorPts[4]=donorPts_[4];
	     //  donorPts[5]=donorPts_[5];

	 // cout << "coucou" <<" "<<irac<<param_intt[NoD][LEVEL]<< param_intt[NoR][LEVEL]<< endl;

	     //cout <<donorPts[0] <<" "<< donorPts[1]<< " "<<  donorPts[2]<<" "<< donorPts[3]<<" "<< donorPts[4]<<" "<< donorPts[5]<< endl;
	      if (nstep%cycle == 1)

		     {

		       //cout << "NoD= " << NoD << endl;

		       // if (ithread==3) 
		       // {

			   //  cout << donorPts[0] <<" "<< donorPts[1] <<" "<< donorPts[2] <<" "<<donorPts[3] << " "<< irac<< endl; 

		       // }

		       //#pragma omp barrier


		       /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn 
		       interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre,iptstk + irac*taille + 1*5*taillefenetre, iptro_p1[NoD],dir,taillefenetre,nstep,NoD);
		     
		     }
	     

		   if (nstep%cycle == cycle/2 + cycle/4) 
		     {

		       //cout << "coucou cell fictives" << "  "<<cycle <<"  "<< NoD<<endl;
		       //#pragma omp barrier

		       //E_Int ind=2;
		       //copy_rk3localpara_(param_intt[NoD], donorPts , donorPts_, iptro[NoD], iptstk + irac*taille + 0*5*taillefenetre,ind,taillefenetre);

		       
#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

		       E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 ];
		       E_Int pos;
		       pos  = ipt_param_int[ shift_rac + nrac*6  ]; E_Int* donorPts2   = ipt_param_int +  pos;
		       pos  = ipt_param_int[ shift_rac + nrac*12 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor

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
		       	     { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); pt_fin=pt_fin-1; }
		       	   else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk;  pt_fin=pt_fin-1; } 

		       	   //if (ithread==8)
		       	   // {
		       	   //   cout << pt_deb <<" "<< pt_fin<<" "<<nbRcvPts<< endl;
		       	   // }
		       	 }
		       else
		       	 {
			   pt_deb=0;
			   pt_fin=nbRcvPts-1;
			   nbRcvPts = pt_fin - pt_deb +1;
			   //cout << "coucou" << endl;
			   }

		       //cout << irac << endl;

		       remp_cellfictivespara_(param_intt[NoD],param_intt[NoR], iptro[NoD], iptro[NoR],donorPts2,rcvPts,nbRcvPts,pt_deb,pt_fin);



		     }


      }

     else if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 
       {

	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	     E_Int donorPts_[6]; 
             E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
             donorPts_[0] = ipt_param_int[adr + 7];
             donorPts_[1] = ipt_param_int[adr + 8] ;
             donorPts_[2] = ipt_param_int[adr + 9];
             donorPts_[3] = ipt_param_int[adr +10];
             donorPts_[4] = ipt_param_int[adr +11];
             donorPts_[5] = ipt_param_int[adr +12];
             E_Int dir    = ipt_param_int[adr +13];

	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

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

	     //  donorPts[0]=donorPts_[0];
	     // donorPts[1]=donorPts_[1];
	     //  donorPts[2]=donorPts_[2];
	     //  donorPts[3]=donorPts_[3];
	     //  donorPts[4]=donorPts_[4];
	     // donorPts[5]=donorPts_[5];


	 //cout << "coucou" <<" "<<irac<<param_intt[NoD][LEVEL]<< param_intt[NoR][LEVEL]<< endl;
	 if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
	   //if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1)

	   {
	     //cout << "coucou cell fictives" << endl;

	     //#pragma omp barrier

	     /// Recuperation de yinterpolé en position 0 dans le tableau de stockage du raccord
	     //E_Int ind=2;    
	     //copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_ , iptro[NoD], iptstk + irac*taille + 0*taillefenetre,ind,taillefenetre); 

	     /// Ecriture de yinterpolé dans rop
	     E_Float coeff=2.0 ;//cout << NoD << endl;		        	       
	     interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, iptstk + irac*taille + 0*taillefenetre, iptdrodmstk + irac*taille + 0*taillefenetre,iptro[NoD],taillefenetre,coeff);  

#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

	     E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];
	     
	     // //pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
	     // //pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
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
	     	   else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk;  pt_fin=pt_fin-1; } 
	       }
	       else
	     	 {
		   pt_deb=0;
		   pt_fin=nbRcvPts-1;
		   nbRcvPts = pt_fin - pt_deb +1;
		 }

	     
		    //cout << irac << endl;
	      remp_cellfictivespara_(param_intt[NoD],param_intt[NoR], iptro[NoD], iptro[NoR],donorPts2,rcvPts,nbRcvPts,pt_deb,pt_fin);


	   }
	     
       }

	      

    }

  } // fin zone omp






  /* 
 

  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

    {

      E_Int shift_rac =  ech + 2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];

       
     E_Int cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];

     if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse 
       {

	 // cout << "coucou" <<" "<<irac<<param_intt[NoD][LEVEL]<< param_intt[NoR][LEVEL]<< endl;
	      if (nstep%cycle == 1)

		     {

		       E_Int donorPts[6];  E_Int idir = 2;
		       donorPts[0] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 7];
		       donorPts[1] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 8] ;
		       donorPts[2] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 9];
		       donorPts[3] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 10];
		       donorPts[4] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 11];
		       donorPts[5] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 12];
		       int dir = ipt_param_int[ech + 2 + nrac*16 + 14*irac + 13];

		       donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

	       
		       E_Int taillefenetre;
		       taillefenetre = (donorPts[1] - donorPts[0] + 1)*(donorPts[3] - donorPts[2] + 1)*(donorPts[5] - donorPts[4] + 1);

		       /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn 
		       interp_rk3local3_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts, iptstk + irac*taille + 0*taillefenetre,iptstk + irac*taille + 1*5*taillefenetre, iptro_p1[NoD],dir,taillefenetre,nstep);
		     
		     }



		   if (nstep%cycle == cycle/2 + cycle/4) 
		     {
		       /// Interpolation de y6
		       E_Int donorPts[6];  E_Int idir = 2;
		       donorPts[0] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 7];
		       donorPts[1] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 8];
		       donorPts[2] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 9];
		       donorPts[3] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 10];
		       donorPts[4] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 11];
		       donorPts[5] =  ipt_param_int[ech + 2 + nrac*16 + 14*irac + 12];
		       int dir = ipt_param_int[ech + 2 + nrac*16 + 14*irac + 13];

		       donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;

		       E_Int taillefenetre;
		       taillefenetre = (donorPts[1] - donorPts[0] + 1)*(donorPts[3] - donorPts[2] + 1)*(donorPts[5] - donorPts[4] + 1);

		       E_Int ind=2;
		       copy_rk3local_(param_intt[NoD], donorPts , iptro[NoD], iptstk + irac*taille + 0*5*taillefenetre,ind,taillefenetre);

		       E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];
		       E_Int pos;
		       //pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
		       //pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
		       pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts2   = ipt_param_int +  pos;
		       pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
		       //pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

		       //cout << irac << endl;

		       remp_cellfictives_(param_intt[NoD],param_intt[NoR], iptro[NoD], iptro[NoR],donorPts2,rcvPts,nbRcvPts);

		     }



       }


     else if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 
       {

	 //cout << "coucou" <<" "<<irac<<param_intt[NoD][LEVEL]<< param_intt[NoR][LEVEL]<< endl;
	 if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
	   //if (nstep%cycle==cycle-2 and (nstep/cycle)%2==1)

	   {


	     E_Int pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; 
	     E_Int donorPts[6]; 
	     donorPts[0] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 7];
             donorPts[1] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 8] ;
             donorPts[2] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 9];
             donorPts[3] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 10];
             donorPts[4] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 11];
             donorPts[5] =  ipt_param_int[ech + 2 + nrac*14 + 14*irac + 12];
	     int dir = ipt_param_int[ech + 2 + nrac*14 + 14*irac + 13];


	     donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*1;


	     //cout << "irac= " <<irac<< endl;
	     E_Int taillefenetre;
	     taillefenetre = (donorPts[1] - donorPts[0] + 1)*(donorPts[3] - donorPts[2] + 1)*(donorPts[5] - donorPts[4] + 1); 


	     /// Recuperation de yinterpolé en position 0 dans le tableau de stockage du raccord
	     E_Int ind=2;    
	     copy_rk3local_(param_intt[NoD], donorPts , iptro[NoD], iptstk + irac*taille + 0*taillefenetre,ind,taillefenetre); 


	     E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];
	     
	     //pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int +  pos;
	     //pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int +  pos;
	     pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts2   = ipt_param_int +  pos;
	     pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int +  pos;   // donor et receveur inverser car storage donor
	     //pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

	     //cout << irac << endl;

	     remp_cellfictives_(param_intt[NoD],param_intt[NoR], iptro[NoD], iptro[NoR],donorPts2,rcvPts,nbRcvPts);

	   }



       }

    }

 
*/
 














   
 

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

  E_Int taille=2000000/nrac;

  for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

    {

      E_Int timelevel = ipt_param_int[ ech +3 ]; 
      //E_Int shift_rac =  ech + 2 + irac;
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5  ]; // Numero zone donneuse du irac concerné
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 ]; // Numero zone receveuse du irac concerné
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 ];

       
     E_Int cycle = param_intt[NoD][NSSITER]/param_intt[NoD][LEVEL];
     E_Int cycleR = param_intt[NoR][NSSITER]/param_intt[NoR][LEVEL];


     if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
       {

	 E_Int pos;E_Int ind;
	     pos  = ipt_param_int[ shift_rac + nrac*6 ]; 
	     E_Int donorPts_[6]; 
             E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
             donorPts_[0] = ipt_param_int[adr + 0];
             donorPts_[1] = ipt_param_int[adr + 1] ;
             donorPts_[2] = ipt_param_int[adr + 2];
             donorPts_[3] = ipt_param_int[adr + 3];
             donorPts_[4] = ipt_param_int[adr + 4];
             donorPts_[5] = ipt_param_int[adr + 5];
             E_Int dir    = ipt_param_int[adr + 6];

	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

	     //cout <<  donorPts_[0]<<" "<< donorPts_[1] <<" "<< donorPts_[2]<<" "<< donorPts_[3] <<" "<< donorPts_[4]<<" "<< donorPts_[5]<<endl;

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
		 //cout << "ind= " << ind << endl;

		 //cout << "coucou_cons" << endl;
		 conservrk3local3para_(param_intt[NoR],donorPts,donorPts_,iptdrodm + shift_zone[NoR],iptcoe+shift_coe[NoR], iptcstk+irac*taille,taillefenetre,nstep,ind);

	       }


	     pos  = ipt_param_int[ shift_rac + nrac*6  ]; 
	      
             E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
             donorPts_[0] = ipt_param_int[adr + 7];
             donorPts_[1] = ipt_param_int[adr + 8] ;
             donorPts_[2] = ipt_param_int[adr + 9];
             donorPts_[3] = ipt_param_int[adr +10];
             donorPts_[4] = ipt_param_int[adr +11];
             donorPts_[5] = ipt_param_int[adr +12];
             dir          = ipt_param_int[adr +13];

	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];
	     //cout << "Pts receveurs : "<<  donorPts[0] <<" "<< donorPts[1] <<" "<< donorPts[2] <<" "<< donorPts[3] << endl;	     
	     taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 

	     
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
		conservrk3local3para_(param_intt[NoD],donorPts,donorPts_,iptdrodm + shift_zone[NoD],iptcoe+shift_coe[NoD], iptcstk+irac*taille+5*taillefenetre,taillefenetre,nstep,ind);

	       }

             E_Int adr = ech + 4 + timelevel*2 + nrac*18 + 14*irac;
             donorPts_[0] = ipt_param_int[adr + 7];
             donorPts_[1] = ipt_param_int[adr + 8] ;
             donorPts_[2] = ipt_param_int[adr + 9];
             donorPts_[3] = ipt_param_int[adr +10];
             donorPts_[4] = ipt_param_int[adr +11];
             donorPts_[5] = ipt_param_int[adr +12];
             dir          = ipt_param_int[adr +13];

	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
	     donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

	     taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 


	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                              donorPts_,
                              topology, ipt_ind_dm_omp_thread);

	     donorPts[0]=ipt_ind_dm_omp_thread[0];
	     donorPts[1]=ipt_ind_dm_omp_thread[1];
	     donorPts[2]=ipt_ind_dm_omp_thread[2];
	     donorPts[3]=ipt_ind_dm_omp_thread[3];
	     donorPts[4]=ipt_ind_dm_omp_thread[4];
	     donorPts[5]=ipt_ind_dm_omp_thread[5];

	     E_Int* transfo;
	     transfo = &ipt_param_int[adr + 14];
	       
             E_Int pt1 =  ipt_param_int[adr + 17];
             E_Int pt2 =  ipt_param_int[adr + 18];
             E_Int pt3 =  ipt_param_int[adr + 19];



	     if (nstep%cycle==cycle-1) /// Opération pour assurer la conservativité
	       {


#pragma omp barrier // On attend que tous les threads aient calculé leur bilan de flux
		 //cout << "zone "<< NoD<< endl;
		 //	 conservrk3local4para_(param_intt[NoD],param_intt[NoR],iptcoe+shift_coe[NoD],iptcoe+shift_coe[NoR],param_realt[NoD],donorPts,donorPts_,iptro_p1[NoD],iptcstk+irac*taille+5*taillefenetre,iptcstk+irac*taille,taillefenetre,nstep,dir,pt1,pt2,pt3,transfo);

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
