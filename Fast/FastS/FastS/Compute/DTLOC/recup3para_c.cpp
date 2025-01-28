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

using namespace std;
using namespace K_FLD;


//=============================================================================
// Idem: in place + from zone + tc compact au niveau base 
//=============================================================================
void K_FASTS::recup3para_c(E_Float**& iptro, E_Int*& param_int_tc,
                           E_Int**& param_int, 
                           E_Float*& iptstk, E_Int& nstep, E_Int& taille_tabs, E_Int& nidom)

{

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

	E_Int ech       = param_int_tc[ NoTransfert +shift_graph];
	E_Int nrac      = param_int_tc[ ech +1 ];
        E_Int timelevel = param_int_tc[ ech +3 ];

	for  (E_Int irac=0; irac< nrac; irac++) // Boucle sur les différents raccords (frontières)
	  {
	    E_Int shift_rac = ech + 4 + timelevel*2 + irac;
	    E_Int debut_rac = ech + 4 + timelevel*2 + nrac*18 + 27*irac;
	    E_Int NoD       =  param_int_tc[ shift_rac + nrac*5     ]; // Numero zone donneuse du irac concerné

	    E_Int cycle = param_int[NoD][NSSITER]/param_int_tc[debut_rac +25];

	    if (param_int_tc[debut_rac + 25] > param_int_tc[debut_rac + 24]) /// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse 
	      {
		if (nstep%cycle==cycle/2 and (nstep/cycle)%2==1) 
		  {
#                   include "fenetre_raccord.h"

		    //// Recuperation de y6 en position 1 dans le tableau de stockage du raccord
		    E_Int ind=2;
		    copy_rk3localpara_(param_int[NoD], donorPts , donorPts_, iptro[NoD], iptstk + pos_tab + 1*5*taillefenetre,ind,taillefenetre); 
		  }
	      }
	    else if (param_int_tc[debut_rac + 25] < param_int_tc[debut_rac + 24]) /// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	      {
		if (nstep%cycle==cycle/2 + cycle/4) 	   
		  {
#                   include "fenetre_raccord.h"

		    E_Int ind=2;
		    copy_rk3localpara_(param_int[NoD], donorPts ,donorPts_,  iptro[NoD], iptstk + pos_tab + 2*5*taillefenetre,ind,taillefenetre);
		  }
	      }
	  }

      } //boucle NoTransfert

  } // fin zone omp

  //#pragma omp barrier



}
