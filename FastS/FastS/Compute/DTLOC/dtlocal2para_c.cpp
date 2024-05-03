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
void K_FASTS::dtlocal2para_c(E_Float**& iptro, E_Float**& iptro_p1,
                             E_Int*& param_int_tc, E_Float*& param_real_tc,
                             E_Int**& param_int  , E_Float**& param_real  ,
			     E_Float*& iptdrodm,  E_Float*& iptcoe,  
			     E_Float*& iptstk, E_Float*& iptdrodmstk, E_Float*& iptcstk, E_Int& nstep, E_Int& omp_mode, E_Int& taille_tabs, E_Int& nidom)
{


#ifdef NB_SOCKET
  E_Int nb_socket=NB_SOCKET;
#else
  E_Int nb_socket=1;
#endif
#ifdef CORE_PER_SOCK
  E_Int core_per_socket=CORE_PER_SOCK;
#else
  E_Int core_per_socket=__NUMTHREADS__;
#endif

  if( __NUMTHREADS__ <= core_per_socket)  nb_socket=1;

  //printf(" socket core %d %d \n", nb_socket, core_per_socket);
  
  FldArrayI tab_activ_core_per_socket(nb_socket);   E_Int* activ_core_per_socket = tab_activ_core_per_socket.begin();

  E_Int count = __NUMTHREADS__;
  for (E_Int s = 0; s < nb_socket; s++)
     {   if (count >= core_per_socket) {activ_core_per_socket[s]=  core_per_socket; count -= core_per_socket;}
         else activ_core_per_socket[s] = count;
     }


  E_Int a=0;
  E_Int b=0;


  E_Int shift_zone[nidom];
  E_Int shift_coe[nidom];
  E_Int process = 1;
  
  
  E_Int omp_mode_loc = 0;


  for (E_Int nd = 0; nd < nidom; nd++)
    {
      shift_zone[nd]=a;
      a=a+param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];	 
    }
  for (E_Int nd = 0; nd < nidom; nd++)
    {
      shift_coe[nd]=b;
      b=b+param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];	 
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
    E_Int Nbre_thread_actif_loc = Nbre_thread_actif;
    E_Int ithread_loc = ithread;

    E_Int topology[3];
    topology[0]=0;
    topology[1]=0;
    topology[2]=0;


    for (E_Int NoTransfert = 1; NoTransfert < param_int_tc[0] + 1; NoTransfert++)

      {

    	E_Int nbcomIBC = param_int_tc[1];
    	E_Int nbcomID  = param_int_tc[2+nbcomIBC];
  
    	E_Int shift_graph = nbcomIBC + nbcomID + 2;

    	E_Int ech            = param_int_tc[ NoTransfert +shift_graph];
 
    	E_Int nrac = param_int_tc[ ech +1 ];

	for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

	  {

	    E_Int timelevel = param_int_tc[ ech +3 ]; 
	    E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

	    E_Int NoD      =  param_int_tc[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	    E_Int NoR      =  param_int_tc[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
	    E_Int nvars_loc=  param_int_tc[ shift_rac + nrac*13 +1 ];
 
     
            //Si omp_mode=1, on modifie la distribution du travail pour ameliorer le numa
	    //#           include "distrib_omp.h"

	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;

	    E_Int ibcType = param_int_tc[shift_rac + nrac * 3];

    
	    E_Int cycle = param_int[NoD][NSSITER]/param_int_tc[debut_rac + 25];

	    
	    if (param_int_tc[debut_rac + 25] > param_int_tc[debut_rac + 24])  /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	      {

                //cout << param_int_tc[debut_rac + 25] << "  " <<  param_int_tc[debut_rac + 24] << endl;
												 
		if (nstep%cycle==1 and (nstep/cycle)%2==0) /// La 1ere sous-iteration du cycle
		  {

		    
#                   include "fenetre_raccord.h"
           	    //// stockage des flux en vue de l'interpolation pour la zone de + gd pas de temps (ici NoR) en position 0 dans le tableau de stockage des flux (raccord)
		    E_Int ind=1;
		    copyflux_rk3localpara_(param_int[NoD],donorPts, donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + pos_tab + 0*taillefenetre, ind, taillefenetre);

		    
		    /// stockage de yn en position 0 dans le tableau de stockage du raccord
		    copy_rk3localpara_(param_int[NoD], donorPts,donorPts_, iptro[NoD], iptstk + pos_tab + 0*taillefenetre,ind,taillefenetre);


		  }

		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==0 ) /// Recuperation des valeurs stockées en position 0 (yn ou yinterpolé suivant la parité du cycle) et stockage de y2 ou y6
		  {   
#                   include "fenetre_raccord.h"
		    //// Recuperation et transformation de f(yn) + coeff*f(y1) pour obtenir alpha*f(yn) + beta*f(y1)
		    E_Int ind=2;
		    copyflux_rk3localpara_(param_int[NoD], donorPts , donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + pos_tab + 0*taillefenetre, ind, taillefenetre);
		  }

		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
		  {
#                   include "fenetre_raccord.h"
		    /// Stockage de y6  en position 1 dans le tableau de stockage du raccord
		    E_Int ind=1;
		    copy_rk3localpara_(param_int[NoD], donorPts ,donorPts_, iptro[NoD], iptstk + pos_tab + 1*5*taillefenetre,ind,taillefenetre); 
		  }


	      }

	    else if (param_int_tc[debut_rac + 25] < param_int_tc[debut_rac + 24]) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	      {

		if (nstep%cycle == 1)
		  {

		
#                   include "fenetre_raccord.h"

		    /// Stockage de yn en position 0 dans le tableau de stocakge du raccord 
		    E_Int ind=1;
		    copy_rk3localpara_(param_int[NoD], donorPts , donorPts_, iptro[NoD], iptstk + pos_tab + 0*taillefenetre,ind,taillefenetre); 


		    /// Stockage de y3 en position 1 dans le tableau de stockage du raccord 
		    copy_rk3localpara_(param_int[NoD], donorPts , donorPts_, iptro_p1[NoD], iptstk + pos_tab + 1*5*taillefenetre,ind,taillefenetre); 

    
		    /// stockage de drodm (yn)
		    ind=1;
		    copyflux_rk3local2para_(param_int[NoD], donorPts , donorPts_ ,iptdrodm + shift_zone[NoD], iptdrodmstk + pos_tab + 0*taillefenetre,ind,taillefenetre);
		  }
	      
	      
		if (nstep%cycle == cycle/4)
		  {		  

#                   include "fenetre_raccord.h"
		    /// Interpolation de y2
		    interp_rk3local3para_(param_int[NoD], param_real[NoD], iptcoe+shift_coe[NoD], donorPts, donorPts_,
                                          iptstk + pos_tab + 0*taillefenetre, iptstk + pos_tab + 1*5*taillefenetre, iptro[NoD], dir, taillefenetre, nstep, process);

		  }
		if (nstep%cycle == cycle/2-1)
		  {
#                   include "fenetre_raccord.h"
		    /// Recuperation de y3 en position 1 dans le tableau de stockage du raccord 
		    E_Int ind=2;
		    copy_rk3localpara_(param_int[NoD], donorPts ,donorPts_, iptro_p1[NoD], iptstk + pos_tab + 1*5*taillefenetre,ind,taillefenetre); 

		  }


		if (nstep%cycle == cycle/2)

		  {
#                   include "fenetre_raccord.h"
		    /// Switch pointeurs

		    /// Stockage de y4 en position 4 dans le tableau de stockage du raccord 
		    E_Int ind=1;
		    copy_rk3localpara_(param_int[NoD], donorPts ,donorPts_, iptro[NoD], iptstk + pos_tab + 2*5*taillefenetre,ind,taillefenetre);

		    /// stockage de drodm(y1) pour obtenir alpha*drodm(yn) + beta*drodm(y1)
		    ind=2;
		    copyflux_rk3local2para_(param_int[NoD], donorPts , donorPts_, iptdrodm + shift_zone[NoD], iptdrodmstk + pos_tab + 0*taillefenetre,ind,taillefenetre);


		  }

		if (nstep%cycle == cycle/2 + 1)
		  {
#                   include "fenetre_raccord.h"
		    /// Switch pointeurs

		    //Interpolation de y6
		    interp_rk3local3para_(param_int[NoD], param_real[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_,
                                          iptstk + pos_tab + 0*taillefenetre, iptstk + pos_tab + 1*5*taillefenetre, iptro_p1[NoD], dir, taillefenetre, nstep, process);
		  }

	      }
 

	  }// boucle raccords

      } // boucle NoTransfert

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
    E_Int Nbre_thread_actif_loc = Nbre_thread_actif;
    E_Int ithread_loc = ithread;

    E_Int topology[3];
    topology[0]=0;
    topology[1]=0;
    topology[2]=0;

    //cout << "Nb_transferts= " <<  param_int_tc[0] + 1 << endl; 

    for (E_Int NoTransfert = 1; NoTransfert < param_int_tc[0] + 1; NoTransfert++)

      {

    	E_Int nbcomIBC = param_int_tc[1];
    	E_Int nbcomID  = param_int_tc[2+nbcomIBC];
  
    	E_Int shift_graph = nbcomIBC + nbcomID + 2;

    	//E_Int threadmax_sdm  = __NUMTHREADS__;
    	E_Int ech            = param_int_tc[ NoTransfert +shift_graph];

 
    	E_Int nrac = param_int_tc[ ech +1 ];

	//if (NoTransfert==1){cout << "process, nrac= " << process <<" "<< nrac << endl;}

	for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

	  {
	    //E_Int shift_rac =  ech + 2 + irac;
	    E_Int timelevel = param_int_tc[ ech +3 ]; 
	    E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

	    E_Int NoD      =  param_int_tc[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	    E_Int NoR      =  param_int_tc[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
	    E_Int nvars_loc=  param_int_tc[ shift_rac + nrac*13 +1 ];


	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;

	    E_Int ibcType = param_int_tc[shift_rac + nrac * 3];

	    E_Int cycle;
	    cycle = param_int[NoD][NSSITER]/param_int_tc[debut_rac + 25];




	    if (param_int_tc[debut_rac + 25] < param_int_tc[debut_rac + 24]) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse

	      {




		if (nstep%cycle == 1 )//and irac == 5 and process == 0)
		  {
                    //Si omp_mode=1, on modifie la distribution du travail pour ameliorer le numa
		    //#                   include "distrib_omp.h"
#                   include "fenetre_raccord.h"

		    interp_rk3local3para_(param_int[NoD], param_real[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_, 
                                          iptstk + pos_tab + 0*taillefenetre, iptstk + pos_tab + 1*5*taillefenetre, iptro_p1[NoD], dir, taillefenetre, nstep, process);
		  }
	     
		if (nstep%cycle == cycle/2 + cycle/4)
		  {

		    //#                   include "distrib_omp.h"
#                   include "fenetre_raccord.h"

		    //Interpolation de y6 ds rop
		    E_Float coeff=1.0;
		    interp_rk3local4para_(param_int[NoD], param_real[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_,
                                          iptstk + pos_tab + 0*taillefenetre, iptdrodmstk + pos_tab + 0*taillefenetre, iptro[NoD], taillefenetre, coeff);
		       
#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

		  }


	      }

	    else if (param_int_tc[debut_rac + 25] > param_int_tc[debut_rac + 24]) /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	      {


		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1)
		  {

		    //#                   include "distrib_omp.h"
#                   include "fenetre_raccord.h"
		    /// Ecriture de yinterpolé dans rop
		    E_Float coeff=2.0 ;		        	       
		    interp_rk3local4para_(param_int[NoD], param_real[NoD], iptcoe+shift_coe[NoD], donorPts,donorPts_,
                                          iptstk + pos_tab + 0*taillefenetre, iptdrodmstk + pos_tab + 0*taillefenetre, iptro[NoD], taillefenetre, coeff); 

#pragma omp barrier // On attend que tous les threads aient écrit la solution dans iptro[NoD]

		  }
	      }
	  }
	 } //boucle NoTransfert
  } // fin zone omp



    
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////// Ajout de la conservativite ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


if(param_int[0][ EXPLOCTYPE ] == 1)

  {


    E_Int taille;

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


    for (E_Int NoTransfert = 1; NoTransfert < param_int_tc[0] + 1; NoTransfert++)

      {

    	E_Int nbcomIBC = param_int_tc[1];
    	E_Int nbcomID  = param_int_tc[2+nbcomIBC];
  
    	E_Int shift_graph = nbcomIBC + nbcomID + 2;

    	E_Int ech            = param_int_tc[ NoTransfert +shift_graph];
 
    	E_Int nrac = param_int_tc[ ech +1 ];

	taille=taille_tabs/nrac;


 
	for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)

	  {

	    E_Int timelevel = param_int_tc[ ech +3 ]; 
	    //E_Int shift_rac =  ech + 2 + irac;
	    E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

	    E_Int NoD      =  param_int_tc[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné
	    E_Int NoR      =  param_int_tc[ shift_rac + nrac*11 +1 ]; // Numero zone receveuse du irac concerné
	    E_Int nvars_loc=  param_int_tc[ shift_rac + nrac*13 +1 ];

	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;       

	    E_Int cycle  = param_int[NoD][NSSITER]/param_int_tc[debut_rac + 25];
	    E_Int cycleR = param_int[NoR][NSSITER]/param_int_tc[debut_rac + 24];

 

	    if (param_int_tc[debut_rac + 25] < param_int_tc[debut_rac + 24]) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	      {


		E_Int pos;E_Int ind;
		pos  = param_int_tc[ shift_rac + nrac*6      ]; 
		E_Int donorPts_[6]; 
		donorPts_[0] =  param_int_tc[debut_rac + 0];
		donorPts_[1] =  param_int_tc[debut_rac + 1];
		donorPts_[2] =  param_int_tc[debut_rac + 2];
		donorPts_[3] =  param_int_tc[debut_rac + 3];
		donorPts_[4] =  param_int_tc[debut_rac + 4];
		donorPts_[5] =  param_int_tc[debut_rac + 5];
		int dir = param_int_tc[debut_rac + 6];

		//cout << "dir= " << dir << endl;

		donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
		donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

		//cout << "receveurs_cons= " << donorPts_[0]<<" "<< donorPts_[1] <<" "<< donorPts_[2]<<" "<< donorPts_[3] << endl;



		E_Int taillefenetreR;
		taillefenetreR = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 

		E_Int ipt_ind_dm_omp_thread[6]; 

		indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_int[NoD][ ITYPCP ],
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



		    conservrk3local3para_(param_int[NoR],donorPts,donorPts_,iptdrodm + shift_zone[NoR],iptcoe+shift_coe[NoR], iptcstk+irac*taille,taillefenetreR,nstep,ind);

		  }


		pos  = param_int_tc[ shift_rac + nrac*6      ]; 
	      


		donorPts_[0] =  param_int_tc[debut_rac + 7];
		donorPts_[1] =  param_int_tc[debut_rac + 8] ;
		donorPts_[2] =  param_int_tc[debut_rac + 9];
		donorPts_[3] =  param_int_tc[debut_rac + 10];
		donorPts_[4] =  param_int_tc[debut_rac + 11];
		donorPts_[5] =  param_int_tc[debut_rac + 12];
		dir = param_int_tc[debut_rac + 13];

		//cout << "donneurs= "<<  donorPts_[0] <<" "<< donorPts_[1] <<" "<< donorPts_[2] <<" "<< donorPts_[3] <<" "<<dir<< endl;


		//donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_int_tct[NoD][NIJK]-1-4) ;
		//donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_int_tct[NoD][NIJK]-1-4);


	     
		donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
		donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

		//cout << "donneurs_conserv= "<<  donorPts_[0] <<" "<< donorPts_[1] <<" "<< donorPts_[2] <<" "<< donorPts_[3] << endl;

		E_Int taillefenetreD;	     
		taillefenetreD = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 

	     
		indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_int[NoD][ ITYPCP ],
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
		    conservrk3local3para_(param_int[NoD],donorPts,donorPts_,iptdrodm + shift_zone[NoD],iptcoe+shift_coe[NoD], iptcstk+irac*taille+5*taillefenetreR,taillefenetreD,nstep,ind);

		  }


		donorPts_[0] =  param_int_tc[debut_rac + 7];
		donorPts_[1] =  param_int_tc[debut_rac + 8] ;
		donorPts_[2] =  param_int_tc[debut_rac + 9];
		donorPts_[3] =  param_int_tc[debut_rac + 10];
		donorPts_[4] =  param_int_tc[debut_rac + 11];
		donorPts_[5] =  param_int_tc[debut_rac + 12];
		//int dir = param_int_tc[ech + 4 + timelevel*2 + 1 + nrac*16 + 14*irac + 13];

		//donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_int_tct[NoD][NIJK]-4-1) ;
		//donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = 1 + 0.5*(dir+abs(dir))*(param_int_tct[NoD][NIJK]-4-1);

		taillefenetreD = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 


		indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_int[NoD][ ITYPCP ],
				  donorPts_,
				  topology, ipt_ind_dm_omp_thread);

		donorPts[0]=ipt_ind_dm_omp_thread[0];
		donorPts[1]=ipt_ind_dm_omp_thread[1];
		donorPts[2]=ipt_ind_dm_omp_thread[2];
		donorPts[3]=ipt_ind_dm_omp_thread[3];
		donorPts[4]=ipt_ind_dm_omp_thread[4];
		donorPts[5]=ipt_ind_dm_omp_thread[5];

		//dir = param_int_tc[ech + 4 + timelevel*2 + 1 + nrac*16 + 20*irac + 13];
		//donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
		//donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];


		E_Int donorPts__[6]; 
		donorPts__[0] =  param_int_tc[debut_rac + 0];
		donorPts__[1] =  param_int_tc[debut_rac + 1] ;
		donorPts__[2] =  param_int_tc[debut_rac + 2];
		donorPts__[3] =  param_int_tc[debut_rac + 3];
		donorPts__[4] =  param_int_tc[debut_rac + 4];
		donorPts__[5] =  param_int_tc[debut_rac + 5];
		//dir = param_int_tc[ech + 4 +timelevel*2 + 1 + nrac*16 + 20*irac + 6];
		//donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] +(dir/abs(dir))*1;
		//donorPts_[2*abs(dir)-1-abs(dir-abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))];

		E_Int* transfo;
		transfo = &param_int_tc[debut_rac + 14];
  

		E_Int trans1 =  param_int_tc[debut_rac + 14];
		E_Int trans2 =  param_int_tc[debut_rac + 15];
		E_Int trans3 =  param_int_tc[debut_rac + 16];
	       

		E_Int pt1 =  param_int_tc[debut_rac + 17];
		E_Int pt2 =  param_int_tc[debut_rac + 18];
		E_Int pt3 =  param_int_tc[debut_rac + 19];

		E_Int ratio1 =  param_int_tc[debut_rac + 21];
		E_Int ratio2 =  param_int_tc[debut_rac + 22];
		E_Int ratio3 =  param_int_tc[debut_rac + 23];

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

		    conservrk3local4para_(param_int[NoD],param_int[NoR],iptcoe+shift_coe[NoD],iptcoe+shift_coe[NoR],param_real[NoD],donorPts,donorPts_,donorPts__,iptro_p1[NoD],iptcstk+irac*taille+5*taillefenetreR,iptcstk+irac*taille,taillefenetreR,taillefenetreD,nstep,dir,pt1,pt2,pt3,transfo,ratio1,ratio2,ratio3);

		  }

	      }

	  }

      }// fin boucle processus

  }//fin zone omp

  }// test ajout de la conservativite ou pas


}
