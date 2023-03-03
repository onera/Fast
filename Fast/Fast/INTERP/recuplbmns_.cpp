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
# include "fast.h"
# include "../FastC/FastC/fastc.h"
# include "param_solver.h"
# include <string.h>
# include <omp.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
PyObject* K_FAST::recuplbmns_(PyObject* self, PyObject* args)
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
  E_Int**   param_intt  = new E_Int*[nidomR];
  E_Float** param_realt = new E_Float*[nidomR];
  E_Int**   interp_data = new E_Int*[nidomR];
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

       o            = K_PYTREE::getNodeFromName1(zone     , "FlowSolution#Centers");
       PyObject* t  = K_PYTREE::getNodeFromName1( o       , "Density_P1"          );
       iptro_p1[nd] = K_PYTREE::getValueAF(t, hook);

        o            = K_PYTREE::getNodeFromName1(zone    , "FlowSolution#Centers");
       PyObject* t1  = K_PYTREE::getNodeFromName1( o      , "Density"             );
       iptro[nd]     = K_PYTREE::getValueAF(t1, hook);
     }

  /// Recuperation du tableau de stockage des valeurs
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

    E_Int flag_implicit = 0;

    for (E_Int nd=0; nd < nidomR; nd++) // Boucle sur les differents domaines
      {
         //cout << "Domaine = " << nd << endl;

         // A partir du moment ou il y une zone implicite : tout NS implicite
         if (param_intt[nd][ ITYPCP ] == 1) flag_implicit = 1;
         //cout << "Flag implicit =" << flag_implicit << endl;

         if (param_intt[nd][IFLOW] == 4)
           {
              //cout << "Coucou domaine LBM" << endl;

              E_Int nb_rac_zone = interp_data[nd][0];
              //cout << "Nombre de raccords NS LBM pour la zone " << nb_rac_zone << endl;
              E_Int taille_rac_max = interp_data[nd][1];

              for (E_Int irac_lbmns=0; irac_lbmns < nb_rac_zone; irac_lbmns++)
                {
                   //cout << "Raccord NS LBM pour la zone " << irac_lbmns << endl;
                   E_Int shift_gc = 37 + irac_lbmns*6;

                   E_Int donorPts_[6];
                   donorPts_[0] = interp_data[nd][ shift_gc + 1 ];
                   donorPts_[1] = interp_data[nd][ shift_gc + 2 ];
                   donorPts_[2] = interp_data[nd][ shift_gc + 3 ];
                   donorPts_[3] = interp_data[nd][ shift_gc + 4 ];
                   donorPts_[4] = interp_data[nd][ shift_gc + 5 ];
                   donorPts_[5] = interp_data[nd][ shift_gc + 6 ];
                   //cout << donorPts_[0] <<" "<< donorPts_[1] <<" "<< donorPts_[2] <<" "<< donorPts_[3] <<" "<< donorPts_[4] <<" "<< donorPts_[5] << endl;

                   E_Int ipt_ind_dm_omp_thread[6];

 	           indice_boucle_lu_(nd, ithread_loc, Nbre_thread_actif_loc, param_intt[nd][ ITYPCP ],
                                     donorPts_, topology, ipt_ind_dm_omp_thread);

                   E_Int donorPts[6];
                   donorPts[0] = ipt_ind_dm_omp_thread[0];
                   donorPts[1] = ipt_ind_dm_omp_thread[1];
                   donorPts[2] = ipt_ind_dm_omp_thread[2];
                   donorPts[3] = ipt_ind_dm_omp_thread[3];
                   donorPts[4] = ipt_ind_dm_omp_thread[4];
                   donorPts[5] = ipt_ind_dm_omp_thread[5];
                   //cout << donorPts[0] <<" "<< donorPts[1] <<" "<< donorPts[2] <<" "<< donorPts[3] <<" "<< donorPts[4] <<" "<< donorPts[5] << endl;

                   E_Float* tab1;
                   E_Float* tab2;

                   //Cas explicite : dans ce cas on rapatrie a nstep=3
                   if (nstep==3 and flag_implicit==0) // La derniere sous-iteration et juste apres fillghostcells
                     {
                        //cout << " recuplbmns - Je passe dans nstep==3" << endl;

                        tab1 = iptro_p1[nd];
                        tab2 = iptro[nd];

                        fill_ghostcellspara_(param_intt[nd], donorPts, donorPts_, tab1, tab2);

#pragma omp barrier

                     }

                   if (nstep==nitmax-1 and flag_implicit==1)
                     {
                       //cout << "Bonjour normalement rapatrie icite" << endl;

                       tab1 = iptro_p1[nd];
                       tab2 = iptro[nd];

                       fill_ghostcellspara_(param_intt[nd], donorPts, donorPts_, tab1, tab2);

#pragma omp barrier

                     }


                }
           }
      }

	// for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
	//   {
  //
	//     //E_Int shift_rac =  ech + 2 + irac;
        //     E_Int timelevel = ipt_param_int[ ech +3 ];
	//     E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
  //
	//     E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	//     E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
  //
	//     E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
  //
	//     E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
	//     E_Int cycle;
  //
	//     cycle = param_intt[NoD][NSSITER]/ipt_param_int[debut_rac +25];
  //
	//     E_Float* tab1;
	//     E_Float* tab2;
	//     E_Float* tab3;
  //  	    E_Int ind;
  //
	//     //if (param_intt[NoD][LEVEL] > param_intt[NoR][LEVEL])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	//     if (ipt_param_int[debut_rac + 25] > ipt_param_int[debut_rac + 24] and ibcType < 0)
  //
	//       {
	// 	E_Int donorPts_[6];  E_Int idir = 2;
  //
	// 	donorPts_[0] =  ipt_param_int[debut_rac + 7];
	// 	donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
	// 	donorPts_[2] =  ipt_param_int[debut_rac + 9];
	// 	donorPts_[3] =  ipt_param_int[debut_rac + 10];
	// 	donorPts_[4] =  ipt_param_int[debut_rac + 11];
	// 	donorPts_[5] =  ipt_param_int[debut_rac + 12];
	// 	int dir = ipt_param_int[debut_rac + 13];
	// 	E_Int profondeur = ipt_param_int[debut_rac + 20];
  //
	// 	donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);
  //
	// 	E_Int taillefenetre;
	// 	taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);
  //
	// 	E_Int ipt_ind_dm_omp_thread[6];
  //
	// 	indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
	// 			  donorPts_,
	// 			  topology, ipt_ind_dm_omp_thread);
  //
	// 	if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
	// 	  {
	// 	    //// Recuperation de y6 en position 1 dans le tableau de stockage du raccord
	// 	    ind=2;
  //                   tab1 = iptro[NoD];
  //                   tab2 = iptstk + irac*taille_stk + param_intt[NoD][NEQ]*taillefenetre;
  //
	// 	    copy_rk3localpara_(param_intt[NoD], ipt_ind_dm_omp_thread, donorPts_, tab1, tab2 , ind, taillefenetre);
	// 	  }
	//       }
	//     //else if (param_intt[NoD][LEVEL] < param_intt[NoR][LEVEL]) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	//     else if (ipt_param_int[debut_rac + 25] < ipt_param_int[debut_rac + 24] and ibcType < 0 )
	//       {
  //
	// 	E_Int donorPts_[6];  E_Int idir = 2;
  //
	// 	donorPts_[0] =  ipt_param_int[debut_rac + 7];
	// 	donorPts_[1] =  ipt_param_int[debut_rac + 8] ;
	// 	donorPts_[2] =  ipt_param_int[debut_rac + 9];
	// 	donorPts_[3] =  ipt_param_int[debut_rac + 10];
	// 	donorPts_[4] =  ipt_param_int[debut_rac + 11];
	// 	donorPts_[5] =  ipt_param_int[debut_rac + 12];
	// 	int dir = ipt_param_int[debut_rac + 13];
	// 	E_Int profondeur = ipt_param_int[debut_rac + 20];
  //
  //               //extrude la fenetre sur profondeur en fonction de la direction du raccord: codage tres sioux.....
	// 	donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);
  //
	// 	E_Int taillefenetre;
	// 	taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1);
  //
  //
	// 	E_Int ipt_ind_dm_omp_thread[6];
  //
	// 	indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_intt[NoD][ ITYPCP ],
	// 			  donorPts_,
	// 			  topology, ipt_ind_dm_omp_thread);
  //
	// 	if (nstep%cycle==cycle/2 + cycle/4)
	// 	  {
	// 	    ind=2;
  //                   tab1 = iptro[NoD];
  //                   tab2 = iptstk + irac*taille_stk + 2*param_intt[NoD][NEQ]*taillefenetre;
  //
	// 	    copy_rk3localpara_(param_intt[NoD], ipt_ind_dm_omp_thread, donorPts_, tab1, tab2, ind, taillefenetre);
	// 	  }
	//       }
	//   }
  } // fin zone omp


 RELEASESHAREDN( interpArray  , stk );
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
