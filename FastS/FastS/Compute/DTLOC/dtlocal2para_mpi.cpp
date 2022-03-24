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
PyObject* K_FASTS::dtlocal2para_mpi(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *work;
  PyObject *pyParam_int, *pyParam_real;
  E_Int loc, nstep, vartype;
  E_Int omp_mode, NoTransfert, process;

  if (!PYPARSETUPLE(args,
                    "OOOOOlllll", "OOOOOiiiii",
                    "OOOOOlllll", "OOOOOiiiii",
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work,
                    &vartype, &nstep, &omp_mode, &NoTransfert, &process))
  {
      return NULL;
  }


  vector<PyArrayObject*> hook;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);
  

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
  PyObject* dtlocArray = PyDict_GetItemString(work,"tab_dtloc"); FldArrayF* stk;
  K_NUMPY::getFromNumpyArray(dtlocArray, stk, true); E_Float* iptstk = stk->begin();

  E_Int stk_size    = stk[0].getSize();
  E_Int taille_tabs = stk_size/5;

  E_Float* iptdrodmstk = iptstk+ taille_tabs*3;
  E_Float* iptcstk     = iptstk+ taille_tabs*4;
  
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

      E_Int taille_stk       = taille_tabs*3;
      E_Int taille_drodmstk  = taille_tabs;
      E_Int taille_cstck     = taille_tabs;


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



    for (E_Int NoTransfert = 1; NoTransfert < ipt_param_int[0] + 1; NoTransfert++)

      {


    	E_Int nbcomIBC = ipt_param_int[1];
    	E_Int nbcomID  = ipt_param_int[2+nbcomIBC];
  
    	E_Int shift_graph = nbcomIBC + nbcomID + 2;

    	E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];
 
    	E_Int nrac = ipt_param_int[ ech +1 ];

	for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
	  {

	    E_Int timelevel = ipt_param_int[ ech +3 ]; 
	    E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

	    E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	    E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
	    E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];
      
	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
	    E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];

	    E_Int cycle;
	    cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];

            /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	    if (ipt_param_int[debut_rac + 25] > ipt_param_int[debut_rac + 24] )

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

		E_Int pos_tab = ipt_param_int[debut_rac + 26];


		E_Float* tab1;   
		E_Float* tab2;   
		E_Float* tab3;   
		E_Int ind;
		if (nstep%cycle==1 and (nstep/cycle)%2==0) /// La 1ere sous-iteration du cycle
		  {
		    //// stockage des flux en vue de l'interpolation pour la zone de + gd pas de temps (ici NoR) en position 0 dans le tableau de stockage des flux (raccord)
		    ind=1;
		    tab1 = iptdrodm + shift_zone[NoD];
		    tab2 = iptdrodmstk + pos_tab + 0*taillefenetre;
                 
		    copyflux_rk3localpara_(param_intt[NoD],donorPts, donorPts_, tab1, tab2, ind, taillefenetre);

		    /// stockage de yn en position 0 dans le tableau de stockage du raccord
		    tab1 = iptro[NoD];
		    tab2 = iptstk + pos_tab + 0*taillefenetre;

		    copy_rk3localpara_(   param_intt[NoD], donorPts, donorPts_, tab1, tab2, ind, taillefenetre);
		  }
	 
		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==0 ) /// Recuperation des valeurs stockées en position 0 (yn ou yinterpolé suivant la parité du cycle) et stockage de y2 ou y6
		  { 


		    //// Recuperation et transformation de f(yn) + coeff*f(y1) pour obtenir alpha*f(yn) + beta*f(y1)
		    ind=2;
		    tab1 = iptdrodm + shift_zone[NoD];
		    tab2 = iptdrodmstk + pos_tab + 0*taillefenetre;
		    copyflux_rk3localpara_(param_intt[NoD], donorPts , donorPts_, tab1, tab2, ind, taillefenetre);
		  }

		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
		  {
		    /// Stockage de y6  en position 1 dans le tableau de stockage du raccord
		    ind=1;
		    tab1 = iptro[NoD];
		    tab2 = iptstk + pos_tab + param_intt[NoD][NEQ]*taillefenetre;
		    copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, tab1, tab2, ind, taillefenetre);
		  }
	      }

	    else if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24] ) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
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

		E_Int pos_tab = ipt_param_int[debut_rac + 26];
		

		E_Float* tab1;   
		E_Float* tab2;   
		E_Float* tab3;   
		E_Int ind;
		if (nstep%cycle == 1)
		  {

		    //cout << donorPts[0]<<" "<< donorPts[1]<<"  "<< donorPts[2]<<"  "<< donorPts[3]<<"  "<< donorPts[4]<<"  "<< donorPts[5]<<"  "<< NoD << endl; 

		    /// Stockage de yn en position 0 dans le tableau de stocakge du raccord 
		    ind =1;
		    tab1 = iptro[NoD];
		    tab2 = iptstk + pos_tab;

		    copy_rk3localpara_(param_intt[NoD], donorPts , donorPts_, tab1, tab2 , ind, taillefenetre); 

		    /// Stockage de y3 en position 1 dans le tableau de stockage du raccord 
		    tab1 = iptro_p1[NoD];
		    tab2 = iptstk + pos_tab  + 1*param_intt[NoD][NEQ]*taillefenetre;

		    copy_rk3localpara_(param_intt[NoD], donorPts , donorPts_, tab1, tab2 , ind, taillefenetre); 
    
		    /// stockage de drodm (yn)
		    ind=1;
		    tab1 = iptdrodm + shift_zone[NoD];
		    tab2 = iptdrodmstk + pos_tab;

		    copyflux_rk3local2para_(param_intt[NoD], donorPts , donorPts_ , tab1, tab2, ind, taillefenetre);
		  }
	      
	      
		/// Interpolation de y2
		if (nstep%cycle == cycle/4)
		  {		       
		    tab1 = iptstk + pos_tab ;
		    tab2 = iptstk + pos_tab + param_intt[NoD][NEQ]*taillefenetre;
		    tab3 = iptro[NoD];

		    interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts, donorPts_, tab1, tab2, tab3, dir, taillefenetre, nstep, NoD);
		  }
	      
		/// Recuperation de y3 en position 1 dans le tableau de stockage du raccord 
		if (nstep%cycle == cycle/2-1)
		  {
		    ind=2;
		    tab1 =  iptro_p1[NoD];
		    tab2 =  iptstk + pos_tab + param_intt[NoD][NEQ]*taillefenetre;

		    copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, tab1, tab2, ind, taillefenetre); 
		  }

		if (nstep%cycle == cycle/2)
		  {
		    /// Stockage de y4 en position 4 dans le tableau de stockage du raccord 
		    ind=1;
		    tab1 = iptro[NoD];
		    tab2 = iptstk + pos_tab + 2*param_intt[NoD][NEQ]*taillefenetre;

		    copy_rk3localpara_(param_intt[NoD], donorPts ,donorPts_, tab1, tab2, ind, taillefenetre);

		    /// stockage de drodm(y1) pour obtenir alpha*drodm(yn) + beta*drodm(y1)
		    ind=2;
		    tab1 = iptdrodm + shift_zone[NoD];
		    tab2 = iptdrodmstk + pos_tab;
              
		    copyflux_rk3local2para_(param_intt[NoD], donorPts , donorPts_, tab1, tab2, ind, taillefenetre);
		  }

		if (nstep%cycle == cycle/2 + 1)
		  {
		    //Interpolation de y6
		    tab1 = iptstk + pos_tab ;
		    tab2 = iptstk + pos_tab + param_intt[NoD][NEQ]*taillefenetre;
		    tab3 = iptro_p1[NoD];

		    interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts, donorPts_, tab1, tab2, tab3, dir, taillefenetre, nstep, NoD);
		  }
	      }
	  }// boucle raccords
      }// fin zone omp 
}
	 









	


 
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

    //cout << "Nb_transferts= " <<  ipt_param_int[0] + 1 << endl; 

    for (E_Int NoTransfert = 1; NoTransfert < ipt_param_int[0] + 1; NoTransfert++)

      {

    	E_Int nbcomIBC = ipt_param_int[1];
    	E_Int nbcomID  = ipt_param_int[2+nbcomIBC];
  
    	E_Int shift_graph = nbcomIBC + nbcomID + 2;

    	//E_Int threadmax_sdm  = __NUMTHREADS__;
    	E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];

 
    	E_Int nrac = ipt_param_int[ ech +1 ];

	for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

	  {

	    //E_Int shift_rac =  ech + 2 + irac;
	    E_Int timelevel = ipt_param_int[ ech +3 ]; 
	    E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

	    E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	    E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
	    E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];

	    E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];

	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;

	    E_Int cycle;
	    cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];

	    if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24] )/// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
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
		
		E_Int pos_tab = ipt_param_int[debut_rac + 26];

		E_Float* tab1;   
		E_Float* tab2;   
		E_Float* tab3;   
		if (nstep%cycle == 1 )
		  {
		    /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn 
		    tab1 = iptstk + pos_tab ;
		    tab2 = iptstk + pos_tab + param_intt[NoD][NEQ]*taillefenetre;
		    tab3 = iptro_p1[NoD];

		    interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts, donorPts_, tab1, tab2, tab3, dir, taillefenetre, nstep, NoD);
		  }
		if (nstep%cycle == cycle/2 + cycle/4)
		  {
		    //Interpolation de y6 ds rop
		    E_Float coeff=1.0;
		    tab1 = iptstk + pos_tab ;
		    tab2 = iptdrodmstk + pos_tab; 
		    tab3 = iptro[NoD];

		    interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, tab1, tab2, iptro[NoD], taillefenetre, coeff);

		    //// are you sure???
		    //// are you sure???
		    //// are you sure???
		    //// are you sure???
		    //// are you sure???
		    //// are you sure???
#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

		  }


	      }

	    else if (ipt_param_int[debut_rac + 25] > ipt_param_int[debut_rac + 24])/// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 

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

		E_Int pos_tab = ipt_param_int[debut_rac + 26];

		E_Float* tab1;   
		E_Float* tab2;   
		E_Float* tab3;   
		E_Int ind;
		

		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)

		  {
		    /// Ecriture de yinterpolé dans rop
		    E_Float coeff=2.0 ;		        	       
		    tab1 = iptstk + pos_tab ;
		    tab2 = iptdrodmstk + pos_tab; 
		    tab3 = iptro[NoD];

		    interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, tab1, tab2, tab3, taillefenetre, coeff);  

#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

		  }
	      }
	  }
      } // fin zone omp

}
 //RELEASESHAREDN( constk  , cstk );
 //RELEASESHAREDN( drodmstock  , drodmstk );
 RELEASESHAREDN( dtlocArray  , stk );
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );
 RELEASESHAREDN( coeArray,coe);
 RELEASESHAREDN( drodmArray,drodm);


 Py_INCREF(Py_None);
 return Py_None;
}
