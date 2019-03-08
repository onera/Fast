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
# include <math.h>

using namespace std;
using namespace K_FLD;


//=============================================================================
// Idem: in place + from zone + tc compact au niveau base 
//=============================================================================
PyObject* K_FASTS::decoupe_maillage(PyObject* self, PyObject* args)
{
  PyObject *zones;
  E_Int nb_level;

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &zones , &nb_level)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "Oi", &zones , &nb_level)) return NULL; 
#endif

  vector<PyArrayObject*> hook;

  E_Int nidom   = PyList_Size(zones);

  E_Int**   param_int = new E_Int*[nidom];
  E_Float** cfl = new E_Float*[nidom];

  for (E_Int nd = 0; nd < nidom; nd++)
     {   
       PyObject* zone = PyList_GetItem(zones, nd);
 
       PyObject* numerics = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
       param_int[nd]      = K_PYTREE::getValueAI(o, hook);

       o            = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject* t  = K_PYTREE::getNodeFromName1( o        , "CFL");
       cfl[nd] = K_PYTREE::getValueAF(t, hook);

     }

  E_Int taille_tab=35;
  E_Float cfl_max = 0.0;
  E_Float cfl_loc;
  E_Int indice;

  PyObject* list_infos = PyList_New( 0 );

  for (E_Int nd = 0; nd < nidom; nd++)
  
   {

     E_Int ind_loop[6];
     ind_loop[0]=  0;
     ind_loop[1]=  param_int[nd][IJKV];
     ind_loop[2]=  0;
     ind_loop[3]=  param_int[nd][IJKV+1];
     ind_loop[4]=  0;
     ind_loop[5]=  param_int[nd][IJKV+2];
     ind_loop[5]=  10;
      for (E_Int k = ind_loop[4]; k < ind_loop[5]; ++k)
       {
	 for (E_Int j = ind_loop[2]; j < ind_loop[3]; ++j)
	   {
	     for (E_Int i = ind_loop[0]; i < ind_loop[1]; ++i)
	       {

		 indice = (i + param_int[nd][NIJK+3]) + (j+param_int[nd][NIJK+3])*param_int[nd][NIJK] + (k+param_int[nd][NIJK+4])*param_int[nd][NIJK]*param_int[nd][NIJK+1];
		 //indice = i + j*param_int[nd][NIJK] + k*param_int[nd][NIJK]*param_int[nd][NIJK+1];
		 cfl_loc = cfl[nd][indice];
		 if (cfl_loc > cfl_max){cfl_max = cfl_loc;}
		 

	       }	 

	   }	 

       }


   }   


  cfl_max = ceil(cfl_max);
  cfl_max = float(cfl_max);
  cout << cfl_max << endl;

  E_Int nbbloc_i;
  E_Int nbbloc_j;
  E_Int nbbloc_k;
  E_Int cfl_max_loc;
  E_Int indice_tab; 



  for (E_Int nd = 0; nd <nidom ; nd++)

   {




     E_Int ind_loop[6];
     ind_loop[0]=  0;
     ind_loop[1]=  param_int[nd][IJKV];
     ind_loop[2]=  0;
     ind_loop[3]=  param_int[nd][IJKV+1];
     ind_loop[4]=  0;
     ind_loop[5]=  param_int[nd][IJKV+2];
     ind_loop[5]=  10;
     nbbloc_i = ind_loop[1]/taille_tab;
     if (nbbloc_i*taille_tab < ind_loop[1] ){nbbloc_i = nbbloc_i + 1;}

     nbbloc_j = ind_loop[3]/taille_tab;
     if (nbbloc_j*taille_tab < ind_loop[3] ){nbbloc_j = nbbloc_j + 1;}

     nbbloc_k = ind_loop[5]/taille_tab;
     if (nbbloc_k*taille_tab < ind_loop[5] ){nbbloc_k = nbbloc_k + 1;}

     cout << "nbbloc_i= " << nbbloc_i << endl;
     cout << "nbbloc_j= " << nbbloc_j << endl;
     cout << "nbbloc_k= " << nbbloc_k << endl;
      
     E_Int tableau_level[1000000];
     bool  tableau_bool[1000000];
     E_Int valeur;
     E_Int indice_max = 0;

    for (E_Int kbloc=0; kbloc < nbbloc_k; kbloc++)
       {
	 for (E_Int jbloc=0; jbloc < nbbloc_j; jbloc++)
	   {
	     for (E_Int ibloc=0; ibloc < nbbloc_i; ibloc++)
	       {
		cfl_max_loc=0;
 
		 for (E_Int k = kbloc*taille_tab; k < min((kbloc+1)*taille_tab,ind_loop[5]); k++)
		   {
       		     for (E_Int j = jbloc*taille_tab; j < min((jbloc+1)*taille_tab,ind_loop[3]); j++)
		       {	   
			 for (E_Int i = ibloc*taille_tab; i < min((ibloc+1)*taille_tab,ind_loop[1]); i++)
			 {
			   indice = (i + param_int[nd][NIJK+3]) + (j+param_int[nd][NIJK+3])*param_int[nd][NIJK] + (k+param_int[nd][NIJK+4])*param_int[nd][NIJK]*param_int[nd][NIJK+1];
			   //indice = i + j*param_int[nd][NIJK] + k*param_int[nd][NIJK]*param_int[nd][NIJK+1];
			   cfl_loc = cfl[nd][indice];
			   if (cfl_loc > cfl_max_loc){cfl_max_loc = cfl_loc;}
			   
			 }
		       }
		   }

		   indice_tab = ibloc + jbloc*nbbloc_i + kbloc*nbbloc_i*nbbloc_j;
		   valeur = pow(2,nb_level-1) / pow(2, floor(log(cfl_max / cfl_max_loc)/log(2)) );
		   if (valeur == 0){valeur=1;}
		   tableau_level[indice_tab] = valeur;
		   tableau_bool[indice_tab]  = false;
		   if (indice_tab > indice_max){indice_max = indice_tab;}
		   //cout << tableau_level[indice_tab] << " " << indice_tab << endl;

	       }

	   }
       }

     // Verification des rapports entre les niveaux en temps adajcents  
    for (E_Int i=0 ; i < indice_max + 1 ; i++)
      {
	E_Int k =  i/(nbbloc_i*nbbloc_j);
	E_Int j =  (i - k*nbbloc_i*nbbloc_j)/nbbloc_i;
	E_Int i_=  i - j*nbbloc_i - k*nbbloc_i*nbbloc_j;

	// zone de gauche
	if (i > 0 and i_ > 0)
	  {
	    if (float(tableau_level[i-1])/float(tableau_level[i]) < 0.5)
	      {
		tableau_level[i-1]=tableau_level[i]/2;
		cout << "changement niveau i-1" << endl;
	      }
	  }

	// zone de droite
	if (i < indice_max and i_ < nbbloc_i -1)
	  {
	    if (float(tableau_level[i+1])/float(tableau_level[i]) < 0.5)
	      {
		tableau_level[i+1]=tableau_level[i]/2;
		cout << "changement niveau i+1" << endl;
	      }
	    //else if (float(tableau_level[i+1])/float(tableau_level[i]) > 2.0){tableau_level[i+1]=tableau_level[i]*2;}
	  }


	// zone du haut
	if (i+nbbloc_i < indice_max + 1  and j < nbbloc_j-1)
	  {
	    if (float(tableau_level[i+nbbloc_i])/float(tableau_level[i]) < 0.5)
	      {
		tableau_level[i+nbbloc_i]=tableau_level[i]/2;
		cout << "changement niveau j+1" << endl;
	      }
	    //else if (float(tableau_level[i+nbbloc_i])/float(tableau_level[i]) > 2.0){tableau_level[i+nbbcloc_i]=tableau_level[i]*2;}
	  }

	// zone du bas
	if (i - nbbloc_i > 0 and j > 0)
	  {
	    if (float(tableau_level[i-nbbloc_i])/float(tableau_level[i]) < 0.5)
	      {
		tableau_level[i-nbbloc_i]=tableau_level[i]/2;
		cout << "changement niveau j-1" << endl;
	      }
	    //else if (float(tableau_level[i-nbbloc_i])/float(tableau_level[i]) > 2.0){tableau_level[i-nbbloc_i]=tableau_level[i]*2;}
	  }

	// zone devant
	if (k >0)
	  {
	    if (float(tableau_level[i-nbbloc_j*nbbloc_i])/float(tableau_level[i]) < 0.5)
	      {
		tableau_level[i-nbbloc_j*nbbloc_i]=tableau_level[i]/2;
		cout << "changement niveau k-1" << endl;
	      }
	    //else if (float(tableau_level[i+nbbloc_i])/float(tableau_level[i]) > 2.0){tableau_level[i+nbbcloc_i]=tableau_level[i]*2;}
	  }

	// zone derriere
	if ( k < nbbloc_k-1)
	  {
	    if (float(tableau_level[i+nbbloc_j*nbbloc_i])/float(tableau_level[i]) < 0.5)
	      {
		tableau_level[i+nbbloc_j*nbbloc_i]=tableau_level[i]/2;
		cout << "changement niveau k+1" << endl;
	      }
	    //else if (float(tableau_level[i-nbbloc_i])/float(tableau_level[i]) > 2.0){tableau_level[i-nbbloc_i]=tableau_level[i]*2;}
	  }
      }


    //for (E_Int kbloc=0; kbloc < nbbloc_k; kbloc++)
    // {
    //	 for (E_Int jbloc=0; jbloc < nbbloc_j; jbloc++)
    //	   {
    //	     for (E_Int ibloc=0; ibloc < nbbloc_i; ibloc++)
    //	       {
    //		   indice_tab = ibloc + jbloc*nbbloc_i + kbloc*nbbloc_i*nbbloc_j;
		   //cout << tableau_level[indice_tab] << " " << indice_tab << endl;

    //	       }

    //	   }
    //   }



//cout << indice_max << endl;

    E_Int indice_tab_d;
    E_Int indice_tab_diag;
    E_Int indice_tab_h;
    E_Int max_bloc=0;
    E_Int nb_bloc=0;
    E_Int jmin = 10000000;
    E_Int nb_zone = 0;
    E_Int tableau_jmin[2000];
    E_Int tableau_jmax[2000];
    E_Int compteur = 0;
    E_Int indice_tab_deb = 0;

    //PyObject* infos;
    //infos = new PyObject[7];

    for (E_Int ind=0 ; ind < indice_max + 1 ; ind++)

      {
	if (tableau_bool[ind]==false)
	  {

	  
	    E_Int k =  ind/(nbbloc_i*nbbloc_j);
	    E_Int j =  (ind - k*nbbloc_i*nbbloc_j)/nbbloc_i;
	    E_Int i=  ind - j*nbbloc_i - k*nbbloc_i*nbbloc_j;
	    E_Int i_=0;
	    E_Int j_=0;
	    E_Int nbbloc_max = 0;
	    E_Int imax = 0;
	    E_Int imin = 0;
	    E_Int jmin = 0;
	    E_Int jmax = 0;
	    E_Int nbbloc = 0;
	    E_Int jmin_max=100000000;

	    while ( tableau_level[ind] == tableau_level[ind+i_] and i_ < nbbloc_i-1 and tableau_bool[ind+i_] == false)
	      {
		while ( tableau_level[ind + i_+j_*nbbloc_i] == tableau_level[ind] and j_ < nbbloc_j-1 and tableau_bool[ind+ i_+j_*nbbloc_i] == false)
		  {
		    j_ = j_ + 1;
		  }

		i_ = i_ + 1;

		if (j_ < jmin_max){jmin_max = j_;}
		nbbloc = i_*jmin_max;
		//cout << jmin_max << endl;
		if (nbbloc >= nbbloc_max)
		  {
		    nbbloc_max = nbbloc;
		    imin = i;
		    imax = i_;
		    jmin = j;
		    jmax = jmin_max;

		  }
		j_ = 0;
	      }

	    //cout << i << " "<< i+imax-1 <<" "<<j <<" "<< j+jmax-1 << "  " << ind <<  endl;

	    for (E_Int j__ = j ;j__< j+jmax;j__++)
	     {
	       for (E_Int i__= i ;i__< i+imax;i__++)
	    	  {
		    //cout << tableau_level[i__+j__*nbbloc_i] << " "<<  i__+j__*nbbloc_i      << endl;
	    	    tableau_bool[i__ + j__*nbbloc_i] = true ;
	    	  }	


	    }
	    
	    E_Int ideb = imin*taille_tab + 1;
	    E_Int ifin = min((i+imax)*taille_tab + 1,ind_loop[1]+1);
	    E_Int jdeb = jmin*taille_tab + 1;
	    E_Int jfin = min((j+jmax)*taille_tab + 1,ind_loop[3]+1);
	    E_Int kdeb = 1;
	    E_Int kfin = 201;

	    nb_zone = nb_zone + 1;

	    //cout << ideb << " "<< ifin <<" "<<jdeb <<" "<< jfin << "  " << tableau_level[imin + jmin*nbbloc_i] <<  endl;
	    PyObject* l = Py_BuildValue("[i,i,i,i,i,i,i,i]",ideb,ifin,jdeb,jfin,kdeb,kfin,tableau_level[imin + jmin*nbbloc_i],nd);

	    PyList_Append( list_infos, l );

	    Py_DECREF( l );

	  }




      }


 

     cout << nb_zone << endl;
    //cout <<    tableau_jmin[2000] << " " << tableau_jmin[2000] << endl;
   

   }
  return list_infos;
  
  Py_DECREF( list_infos );
  delete [] param_int;
  delete [] cfl;
  //delete [] tableau_level;
 

}






