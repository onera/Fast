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
# include "../FastS/FastS/fastS.h"
# include "fast.h"
# include "../FastASLBM/FastASLBM/fastLBM.h"
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
  E_Int loc, nstep, nitrun, nitmax;
  E_Int omp_mode, NoTransfert, process;

  if (!PYPARSETUPLE(args,
                    "OOOOOllllll", "OOOOOiiiiii",
                    "OOOOOllllll", "OOOOOiiiiii",
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work, &nitmax, &nitrun, &nstep, &omp_mode,  &NoTransfert, &process))
  {
      return NULL;
  }


  vector<PyArrayObject*> hook;

  //Nombre de zones doneuses et receveuses
  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //// Recuperation du tableau param_int de l'arbre t et des veceurs rop et roptmp
  // Transferts se font sur Density et Density_P1
  E_Int** param_intt    = new E_Int*[nidomR];
  E_Float** param_realt = new E_Float*[nidomR];
  E_Int** interp_data   = new E_Int*[nidomR];
  E_Float** iptro_p1    = new E_Float*[nidomR];
  E_Float** iptro       = new E_Float*[nidomR];

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

                       o = K_PYTREE::getNodeFromName1(zone     , "FlowSolution#Centers");
       PyObject*      t  = K_PYTREE::getNodeFromName1( o       , "Density_P1"          );
       iptro_p1[nd]      = K_PYTREE::getValueAF(t, hook);

                       o = K_PYTREE::getNodeFromName1(zone    , "FlowSolution#Centers");
       PyObject*     t1  = K_PYTREE::getNodeFromName1( o      , "Density"             );
       iptro[nd]         = K_PYTREE::getValueAF(t1, hook);
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

  //E_Int shift_dom = 0;

//#pragma omp parallel default(shared) num_threads(1)//private(cycle)
#pragma omp parallel default(shared)
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

    E_Int shift_dom = 0;

    E_Int flag_implicit = 0;

    for  (E_Int nd=0; nd < nidomR; nd++) // Boucle sur les differents domaines
      {
         //cout << "Domaine = " << nd << endl;

         // A partir du moment ou il y une zone implicite : tout NS implicite
         if (param_intt[nd][ ITYPCP ] == 1) flag_implicit = 1;

         if (param_intt[nd][IFLOW] == 4)
           {
              //cout << "Coucou domaine LBM" << endl;

              E_Int nb_rac_zone = interp_data[nd][0];
              //cout << "Nombre de raccords NS LBM pour la zone " << nb_rac_zone << endl;
              E_Int taille_rac_max = interp_data[nd][1];

              for (E_Int irac_lbmns=0; irac_lbmns < nb_rac_zone; irac_lbmns++)
                {
                   //cout << "Raccord NS LBM pour la zone " << irac_lbmns << endl;
                   E_Int shift_donorPts = irac_lbmns*6 + 1;

                   //Points du raccord au complet
                   E_Int donorPts_[6];
                   donorPts_[0] = interp_data[nd][ shift_donorPts + 1 ];
                   donorPts_[1] = interp_data[nd][ shift_donorPts + 2 ];
                   donorPts_[2] = interp_data[nd][ shift_donorPts + 3 ];
                   donorPts_[3] = interp_data[nd][ shift_donorPts + 4 ];
                   donorPts_[4] = interp_data[nd][ shift_donorPts + 5 ];
                   donorPts_[5] = interp_data[nd][ shift_donorPts + 6 ];
                   //printf("donor_ %d %d %d %d %d %d ith= %d \n", donorPts_[0],donorPts_[1],donorPts_[2],donorPts_[3],donorPts_[4],donorPts_[5], ithread );

                   E_Int taillefenetre;
                   taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);
                   taillefenetre = taillefenetre*5*3;
                   //printf("taille fen %d  ith= %d \n", taillefenetre , ithread );

                   E_Int ipt_ind_dm_omp_thread[6];

 	                 indice_boucle_lu_(nd, ithread_loc, Nbre_thread_actif_loc, param_intt[nd][ ITYPCP ],
                                                           donorPts_, topology, ipt_ind_dm_omp_thread);

                   //Points sur lequels le thread va agir
                   E_Int donorPts[6];
                   donorPts[0] = ipt_ind_dm_omp_thread[0];
                   donorPts[1] = ipt_ind_dm_omp_thread[1];
                   donorPts[2] = ipt_ind_dm_omp_thread[2];
                   donorPts[3] = ipt_ind_dm_omp_thread[3];
                   donorPts[4] = ipt_ind_dm_omp_thread[4];
                   donorPts[5] = ipt_ind_dm_omp_thread[5];
                   //printf("DONOR %d %d %d %d %d %d ith= %d \n", donorPts[0],donorPts[1],donorPts[2],donorPts[3],donorPts[4],donorPts[5], ithread );

                   E_Float* tab1;
                   E_Float* tab2;
                   E_Float* tab3;
                   E_Int ind;
                   E_Int sol;

                   if (nstep==0) /// La 1ere sous-iteration et juste avant compute
                     {

                       /// stockage de la sol u_n en position 0 dans le tableau de stockage du raccord

                       ind = 1; //stockage
                       sol = 0; //stockage en position 0
                       tab1 = iptro[nd]; // La solution u_n du domaine LBM se trouve dans Density
                       tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                       //cout << "shift stk = " << shift_dom + irac_lbmns*taille_rac_max << endl;
                       copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);

                       /* DEBUG UTILTY*/
                       //E_Int nistk = (donorPts_[1]-donorPts_[0]) + 1;
                       //E_Int nistk2 = (donorPts_[3]-donorPts_[2]) + 1;
                       //E_Int nistk3 = (donorPts_[5]-donorPts_[4]) + 1;
                       //cout << nistk << " " << nistk2 << " " << nistk3 << endl;
                       //cout << nistk*nistk2*nistk3 << endl;
                       //E_Int ldeb = donorPts[0]+1-donorPts_[0] + nistk*(donorPts[2]-donorPts_[2]) + nistk*nistk2*(donorPts[4]-donorPts_[4]);
                       //cout << "ldeb = " << ldeb << endl;

#pragma omp barrier

                     }

                   if (nstep==1) /// La 1ere sous-iteration apres compute
                     {

                       //cout << "Je passe nstep==1 interplbmns" << endl;
                       /// stockage de la sol un+1 en position 1 dans le tableau de stockage du raccord

                       ind = 1; //stockage
                       sol = 1; //stockage en position 1
                       tab1 = iptro[nd]; // La solution u_n+1 du domaine LBM se trouve dans Density
                       tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;

                       copy_valuespara_( param_intt[nd], donorPts , donorPts_, tab1, tab2, ind, sol, taillefenetre);

                       // Si NS explicite (RK3) alors on fait les interp
                       if (flag_implicit==0)
                          {
                            //cout << "Hello je passe dans interp" << endl;
                            /// interpolation de la sol pour sous-pas RK3 (1 -> 2)

                            tab1 = iptro_p1[nd]; //Pour les transferts, valeur d'interp dans Density_P1
                            //tab1 = iptro[NoD];
                            //cout << nitrun << endl;

                            interp_rk3para_( param_intt[nd], donorPts, donorPts_, tab1, tab2, nstep, nitrun, taillefenetre);
                          }
                       // Si NS implicite (Gear) on met u_n+1 dans Density_P1
                       else
                          {
                            //cout << "Hello je passe dans implicite" << endl;
                            /// Recuperation de la solution un+1 dans Density_P1 pour transfert
                            ind = 2; //Recuperation
                            sol = 1; //Recupere la sol en position 1
                            tab1 = iptro_p1[nd]; //Pour les transferts, valeur d'interp dans Density_P1
                            tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                            //cout << shift_dom + irac_lbmns*taille_rac_max << endl;

                            copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);
                          }

#pragma omp barrier

                     }

                   if (nstep==2)
                     {

                       /// Pas de stockage de la zone LBM car deja en n+1
                       if (flag_implicit==0)
                          {
                            /// interpolation de la sol pour sous-pas RK3 (2 -> 3)
                            tab1 = iptro[nd]; //Pour les transferts, valeur d'interp dans Density
                            tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                            //cout << nitrun << endl;
                            interp_rk3para_( param_intt[nd], donorPts, donorPts_, tab1, tab2, nstep, nitrun, taillefenetre);
                          }
                       // Si NS implicite (Gear) on met u_n+1 dans Density_P1
                       else
                          {
                            //cout << "Hello je passe dans implicite" << endl;
                            /// Recuperation de la solution un+1 dans Density_P1 pour transfert
                            ind = 2; //Recuperation
                            sol = 1; //Recupere la sol en position 1
                            tab1 = iptro_p1[nd]; //Pour les transferts, valeur d'interp dans Density_P1
                            tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                            //cout << shift_dom + irac_lbmns*taille_rac_max << endl;

                            copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);
                          }

#pragma omp barrier

                     }

                   if (nstep==3)
                     {

                       if (flag_implicit==0)
                          {

                            /// Recuperation de la solution un+1 dans Density_P1 pour transfert
                            ind = 2; //Recuperation
                            sol = 1; //Recupere la sol en position 1
                            tab1 = iptro_p1[nd]; //Pour les transferts, valeur d'interp dans Density_P1
                            tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                            //cout << shift_dom + irac_lbmns*taille_rac_max << endl;

                            copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);

                            /// Recuperation de la solution un+1 dans Density
                            tab1 = iptro[nd];

                            copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);

                            /// Stockage de la sol un-1 (en fait on passe un dans un-1)
                            tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                            E_Int sol_i = 0; //La solution en position 0
                            E_Int sol_d = 2; //passe en position 2 dans stk

                            shiftstk_para_( param_intt[nd], donorPts, donorPts_, tab2, sol_i, sol_d, taillefenetre);

                          }
                       // Si NS implicite (Gear) on met u_n+1 dans Density_P1
                       else
                          {
                            //cout << "Hello je passe dans implicite" << endl;
                            /// Recuperation de la solution un+1 dans Density_P1 pour transfert
                            ind = 2; //Recuperation
                            sol = 1; //Recupere la sol en position 1
                            tab1 = iptro_p1[nd]; //Pour les transferts, valeur d'interp dans Density_P1
                            tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                            //cout << shift_dom + irac_lbmns*taille_rac_max << endl;

                            copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);
                          }

#pragma omp barrier

                     }

                   if(nstep>3 and flag_implicit==1)
                     {
                       /// Recuperation de la solution un+1 dans Density_P1 pour transfert
                       ind = 2; //Recuperation
                       sol = 1; //Recupere la sol en position 1
                       tab1 = iptro_p1[nd]; //Pour les transferts, valeur d'interp dans Density_P1
                       tab2 = iptstk + shift_dom + irac_lbmns*taille_rac_max;
                       //cout << shift_dom + irac_lbmns*taille_rac_max << endl;

                       copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);

#pragma omp barrier

                     }

                } //fin boulce raccords

              shift_dom = shift_dom + nb_rac_zone*taille_rac_max;

         } //fin if LBM

       } // boucle raccords
  } // fin zone omp

 RELEASESHAREDN( interpArray  , stk );
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
