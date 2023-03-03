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
# include <iostream>
# include <string.h>
# include "fastS.h"
# include "fast.h"
# include "fastLBM.h"
# include "param_solver.h"
# include <string.h>
# include <omp.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
PyObject* K_FAST::interplbmns_(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *work;
  PyObject *pyParam_int, *pyParam_real;
  E_Int loc, nstep, nitrun, vartype;
  E_Int omp_mode, NoTransfert, process;

  if (!PYPARSETUPLE(args,
                    "OOOOOllllll", "OOOOOiiiiii",
                    "OOOOOllllll", "OOOOOiiiiii",
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work, &vartype, &nitrun, &nstep, &omp_mode,  &NoTransfert, &process))
  {
      return NULL;
  }


  vector<PyArrayObject*> hook;

  //Nombre de zones doneuses et receveuses
  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //// Recuperation du tableau param_int de l'arbre t et des veceurs rop et roptmp
  // Transferts se font sur Density et Density_P1
  E_Int**   param_intt = new E_Int*[nidomR];
  E_Float** param_realt = new E_Float*[nidomR];
  E_Int**   interp_data = new E_Int*[nidomR];
  E_Float** iptro_p1 = new E_Float*[nidomR];
  E_Float** iptro    = new E_Float*[nidomR];

  // On parcourt les zones R et on recupere les tableaux/vecteurs
  for (E_Int nd = 0; nd < nidomR; nd++)
     {
       PyObject* zone = PyList_GetItem(zonesR, nd);

       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"  );
       param_intt[nd]     = K_PYTREE::getValueAI(o, hook);

	               o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real" );
       param_realt[nd]    = K_PYTREE::getValueAF(o, hook);

                       o  = K_PYTREE::getNodeFromName1(numerics, "Interp_data"    );
       if (o != NULL)
       {  interp_data[nd] = K_PYTREE::getValueAI(o, hook);}

       o            = K_PYTREE::getNodeFromName1(zone     , "FlowSolution#Centers");
       PyObject* t  = K_PYTREE::getNodeFromName1( o       , "Density_P1"          );
       iptro_p1[nd] = K_PYTREE::getValueAF(t, hook);

        o            = K_PYTREE::getNodeFromName1(zone    , "FlowSolution#Centers");
       PyObject* t1  = K_PYTREE::getNodeFromName1( o      , "Density"             );
       iptro[nd]     = K_PYTREE::getValueAF(t1, hook);
     }

 /// Recuperation du tableau de stockage des valeurs pour interpolation
  PyObject* interpArray = PyDict_GetItemString(work,"tab_interp"); FldArrayF* stk;
  K_NUMPY::getFromNumpyArray(interpArray, stk, true); E_Float* iptstk = stk->begin();

  E_Int stk_size    = stk[0].getSize();
  //cout << "taille tab stockage " << stk_size << endl;
  E_Int taille_tabs = stk_size/5;

  /*--------------------------------------*/
  /* Extraction tableau int et real de tc */
  /*--------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();

  E_Int nvars;

  if( vartype <= 3 &&  vartype >= 1) nvars =5;
  else                               nvars =6;
  //cout << vartype << nvars << endl;
  //cout << "ipt_param_int[0]= " << ipt_param_int[0] << endl;

  ///// Les shifts pour les zones /////
  E_Int nbcomIBC = ipt_param_int[1];
  E_Int nbcomID  = ipt_param_int[2+nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 2;

  //E_Int threadmax_sdm  = __NUMTHREADS__;
  E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];


  E_Int nrac = ipt_param_int[ ech +1 ];
  //cout << nrac << endl;

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
        E_Int size_rac =  ipt_param_int[ shift_rac ];              // Taille du raccord
        E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
        E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
        E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];
  
        //E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
        E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
  
        //E_Int cycle;
        //cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];
        //cout << size_rac << endl;
        //cout << irac << "  " << NoD << "  " << NoR << "  " << ibcType << endl;
        //cout << param_intt[NoD][IFLOW] << " " << param_intt[NoR][IFLOW] << endl;
  
        /// La zone donneuse est une zone LBM et la zone receveuse est une zone NS
        if (param_intt[NoD][IFLOW]==4 and param_intt[NoR][IFLOW]!=4 )
  
         {
              //cout << "Bonjour" << endl;
              //cout << irac << "  " << NoD << "  " << NoR << "  " << ibcType << endl;
  
  
              E_Int nb_rac_zone = interp_data[NoD][0];
              //cout << "Nombre de raccord NS LBM pour la zone " << nb_rac_zone << endl;
              E_Int taille_rac_max = interp_data[NoD][1];
 
              // AJOUT D'UNE BOUCLE SUR LES RACCORDS LBM NS
              for (E_Int irac_lbmns=0; irac_lbmns < nb_rac_zone; irac_lbmns++)
              {
                 //cout << "   -> Raccord LBM NS pour zone LBM " << irac_lbmns << endl;
  
                 E_Int shift_donorPts = irac_lbmns*6 + 1;
  
                 E_Int donorPts_[6];
 	         donorPts_[0] =  interp_data[NoD][ shift_donorPts+1 ];
                 donorPts_[1] =  interp_data[NoD][ shift_donorPts+2 ];
                 donorPts_[2] =  interp_data[NoD][ shift_donorPts+3 ];
                 donorPts_[3] =  interp_data[NoD][ shift_donorPts+4 ];
                 donorPts_[4] =  interp_data[NoD][ shift_donorPts+5 ];
                 donorPts_[5] =  interp_data[NoD][ shift_donorPts+6 ];
                 //cout << donorPts_[0] <<" "<< donorPts_[1] <<" "<< donorPts_[2] <<" "<< donorPts_[3] <<" "<< donorPts_[4] <<" "<< donorPts_[5] << endl;
  
                 E_Int taillefenetre;
 	         taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);
                 taillefenetre = taillefenetre*5*3;
                 //cout << taillefenetre << endl;
 	         E_Int ipt_ind_dm_omp_thread[6];
  
                 indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
                                   donorPts_, topology, ipt_ind_dm_omp_thread);
  
 	         E_Int donorPts[6];
 	         donorPts[0]=ipt_ind_dm_omp_thread[0];
 	         donorPts[1]=ipt_ind_dm_omp_thread[1];
 	         donorPts[2]=ipt_ind_dm_omp_thread[2];
 	         donorPts[3]=ipt_ind_dm_omp_thread[3];
 	         donorPts[4]=ipt_ind_dm_omp_thread[4];
 	         donorPts[5]=ipt_ind_dm_omp_thread[5];
                 //cout << donorPts[0] << " " << donorPts[1] << " " << donorPts[2] << " " << donorPts[3] << " " << donorPts[4] << " " << donorPts[5] << endl;
  
                 E_Float* tab1;
 	         E_Float* tab2;
 	         E_Float* tab3;
 	         E_Int ind;
                 E_Int sol;
  
 	         if (nstep==0) /// La 1ere sous-iteration et juste avant compute
 	         {
  
                  //cout << "Je passe dans nstep=0" << endl;
 	 	 /// stockage de la sol un en position 0 dans le tableau de stockage du raccord
 	 	 ind = 1; //stockage
                 sol = 0; //stockage en position 0
                 tab1 = iptro[NoD]; // La solution du domaine LBM se trouve dans Density
                 tab2 = iptstk + irac_lbmns*taille_rac_max;
                 //cout << taille_rac_max << " " << taillefenetre << endl;
 	 	 copy_valuespara_( param_intt[NoD], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);
  
 	         }
  
 	         if (nstep==1) /// La 1ere sous-iteration apres compute
 	         {
  
	 	 //cout << "Je passe dans nstep=1" << endl;
                 //// stockage de la sol un+1 en position 1 dans le tableau de stockage du raccord
 	 	 ind = 1; //stockage
                 sol = 1; //stockage en position 1
                 tab1 = iptro[NoD]; // La solution du domaine LBM se trouve dans Density
                 tab2 = iptstk + irac_lbmns*taille_rac_max;
  
 	 	 copy_valuespara_( param_intt[NoD], donorPts , donorPts_, tab1, tab2, ind, sol, taillefenetre);
  
                 /// interpolation de la sol pour sous-pas RK3 (1 -> 2)
                 tab1 = iptro_p1[NoD]; //Pour les transferts, valeur d'interp dans Density_P1
                 //tab1 = iptro[NoD];
                 //E_Int run = 1;
                 //cout << nitrun << endl;
  
                 interp_rk3para_( param_intt[NoD], donorPts, donorPts_, tab1, tab2, nstep, nitrun, taillefenetre);
  
 	         }
  
 	         if (nstep==2)
 	         {
  
                 //cout << "Je passe dans nstep=2" << endl;
                 /// Pas de stockage
                 /// interpolation de la sol pour sous-pas RK3 (2 -> 3)
                 tab1 = iptro[NoD]; //Pour les transferts, valeur d'interp dans Density
                 tab2 = iptstk + irac_lbmns*taille_rac_max;
                 //E_Int run = 1;
                 //cout << nitrun << endl;
 
                 interp_rk3para_( param_intt[NoD], donorPts, donorPts_, tab1, tab2, nstep, nitrun, taillefenetre);
 
 	         }
  
                 if (nstep==3)
                 {
  
                 /// Recuperation de la solution un+1 dans Density_P1 pour transfert
                 ind = 2;
                 sol = 1; //Recupere la sol en position 1
                 tab1 = iptro_p1[NoD]; //Pour les transferts, valeur d'interp dans Density_P1
                 tab2 = iptstk + irac_lbmns*taille_rac_max;
 
                 copy_valuespara_( param_intt[NoD], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);
 
                 /// Recuperation de la solution un+1 dans Density
                 tab1 = iptro[NoD];
 
                 copy_valuespara_( param_intt[NoD], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);
 
                 /// Stockage de la sol un-1 (en fait on passe un dans un-1)
                 tab2 = iptstk + irac_lbmns*taille_rac_max;
                 E_Int sol_i = 0; //La solution en position 0
                 E_Int sol_d = 2; //passe en position 2 dans stk
  
                 shiftstk_para_( param_intt[NoD], donorPts, donorPts_, tab2, sol_i, sol_d, taillefenetre);
  
                 }
 
             }
  
        }
  
     }// boucle raccords
 }// fin zone omp




// #pragma omp parallel default(shared) //private(cycle)
//   {
// #ifdef _OPENMP
//     E_Int  ithread           = omp_get_thread_num() +1;
//     E_Int  Nbre_thread_actif = omp_get_num_threads();
// #else
//     E_Int ithread = 1;
//     E_Int Nbre_thread_actif = 1;
// #endif
//    E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
//    if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;
//
//    E_Int Nbre_thread_actif_loc, ithread_loc;
//    if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
//    else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}
//
//    E_Int topology[3];
//    topology[0]=0;
//    topology[1]=0;
//    topology[2]=0;
//
//
//    E_Float* tab1;
//    E_Float* tab2;
//    E_Float* tab3;
//
//   for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
//
//     {
//
//       //E_Int shift_rac =  ech + 2 + irac;
//       E_Int timelevel = ipt_param_int[ ech +3 ];
//       E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
//
//       E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
//       E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
//       E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ];
//
//       E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
//
//       E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
//
//       E_Int cycle;
//       cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac + 25];
//
//        if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24] )/// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
//        {
//
// 	     E_Int pos;
// 	     pos  = ipt_param_int[ shift_rac + nrac*6      ];
//
// 	     E_Int donorPts_[6];
// 	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
//              donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
//              donorPts_[2] =  ipt_param_int[debut_rac + 9];
//              donorPts_[3] =  ipt_param_int[debut_rac + 10];
//              donorPts_[4] =  ipt_param_int[debut_rac + 11];
//              donorPts_[5] =  ipt_param_int[debut_rac + 12];
// 	     int dir = ipt_param_int[debut_rac + 13];
// 	     E_Int profondeur = ipt_param_int[debut_rac + 20];
//
// 	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);
//
// 	     //cout << "irac= " <<irac<< endl;
// 	     E_Int taillefenetre;
// 	     taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);
//
// 	     E_Int ipt_ind_dm_omp_thread[6];
//
// 	     indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
//                              donorPts_,
//                              topology, ipt_ind_dm_omp_thread);
//
// 	     E_Int donorPts[6];
// 	     donorPts[0]=ipt_ind_dm_omp_thread[0];
// 	     donorPts[1]=ipt_ind_dm_omp_thread[1];
// 	     donorPts[2]=ipt_ind_dm_omp_thread[2];
// 	     donorPts[3]=ipt_ind_dm_omp_thread[3];
// 	     donorPts[4]=ipt_ind_dm_omp_thread[4];
// 	     donorPts[5]=ipt_ind_dm_omp_thread[5];
//
//
// 	     E_Float* tab1;
// 	     E_Float* tab2;
// 	     E_Float* tab3;
// 	     if (nstep%cycle == 1 )
// 	       {
// 		 /// Interpolation de y1 : y1 = 0.5*y3 + 0.5*yn
// 	         tab1 = iptstk + irac*taille_stk ;
//                  tab2 = iptstk + irac*taille_stk + param_intt[NoD][NEQ]*taillefenetre;
//                  tab3 = iptro_p1[NoD];
//
// 	         interp_rk3local3para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts, donorPts_, tab1, tab2, tab3, dir, taillefenetre, nstep, NoD);
// 	       }
// 	     if (nstep%cycle == cycle/2 + cycle/4)
// 	       {
// 		 //Interpolation de y6 ds rop
// 		 E_Float coeff=1.0;
// 	         tab1 = iptstk + irac*taille_stk ;
//                  tab2 = iptdrodmstk + irac*taille_drodmstk;
//                  tab3 = iptro[NoD];
//
// 		 interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, tab1, tab2, iptro[NoD], taillefenetre, coeff);
//
// //// are you sure???
// //// are you sure???
// //// are you sure???
// //// are you sure???
// //// are you sure???
// //// are you sure???
// #pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]
// //// are you sure???
// //// are you sure???
// //// are you sure???
// //// are you sure???
// //// are you sure???
// //// are you sure???
// 	       }
//
//
//       }
//
//      else if (ipt_param_int[debut_rac + 25] > ipt_param_int[debut_rac + 24])/// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
//        {
//
// 	     E_Int pos;
// 	     pos  = ipt_param_int[ shift_rac + nrac*6      ];
// 	     E_Int donorPts_[6];
// 	     donorPts_[0] =  ipt_param_int[debut_rac + 7];
//              donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
//              donorPts_[2] =  ipt_param_int[debut_rac + 9];
//              donorPts_[3] =  ipt_param_int[debut_rac + 10];
//              donorPts_[4] =  ipt_param_int[debut_rac + 11];
//              donorPts_[5] =  ipt_param_int[debut_rac + 12];
// 	     int dir = ipt_param_int[debut_rac + 13];
// 	     E_Int profondeur = ipt_param_int[debut_rac + 20];
// 	     donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);
//
// 	      E_Int taillefenetre;
// 	      taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);
//
// 	     E_Int ipt_ind_dm_omp_thread[6];
//
// 	      indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
//                                donorPts_,
//                                topology, ipt_ind_dm_omp_thread);
// 	      E_Int donorPts[6];
// 	     donorPts[0]=ipt_ind_dm_omp_thread[0];
// 	     donorPts[1]=ipt_ind_dm_omp_thread[1];
// 	     donorPts[2]=ipt_ind_dm_omp_thread[2];
// 	     donorPts[3]=ipt_ind_dm_omp_thread[3];
// 	     donorPts[4]=ipt_ind_dm_omp_thread[4];
// 	     donorPts[5]=ipt_ind_dm_omp_thread[5];
//
//
// 	     if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
//
// 	       {
// 		 /// Ecriture de yinterpolé dans rop
// 		 E_Float coeff=2.0 ;
// 	         tab1 = iptstk + irac*taille_stk ;
//                  tab2 = iptdrodmstk + irac*taille_drodmstk;
//                  tab3 = iptro[NoD];
//
// 		 interp_rk3local4para_(param_intt[NoD], param_realt[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, tab1, tab2, tab3, taillefenetre, coeff);
//
// #pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]
// 	       }
//        }
//     }
//   } // fin zone omp

 //RELEASESHAREDN( constk  , cstk );
 //RELEASESHAREDN( drodmstock  , drodmstk );
 RELEASESHAREDN( interpArray  , stk );
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );
 //RELEASESHAREDN( coeArray,coe);
 //RELEASESHAREDN( drodmArray,drodm);

 Py_INCREF(Py_None);
 return Py_None;
}
