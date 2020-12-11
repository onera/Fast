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



  E_Int taille_bloc=30;
  E_Float cflmax = 0.0;
  E_Float cflmin = 10000000.0;
  E_Float dtmax;
  E_Float dtmin;
  E_Float cfl_loc;
  E_Int indice;
  E_Int exposant_max;

  PyObject* list_infos = PyList_New( 0 );


 E_Int taille_bloc_i[nidom];
 E_Int taille_bloc_j[nidom];
 E_Int taille_bloc_k[nidom];

  

  for (E_Int nd = 0; nd < nidom; nd++) /// Calcul cflmin et cflmax sur le maillage
  
   {


     E_Int cflref_i;
     E_Int i_grad;
     E_Int i_grad_min =  param_int[nd][IJKV];

     E_Int cflref_j;
     E_Int j_grad;
     E_Int j_grad_min =  param_int[nd][IJKV+1];

     E_Int cflref_k;
     E_Int k_grad;
     E_Int k_grad_min =  param_int[nd][IJKV+2];
     

    
     E_Int ind_loop[6];
     ind_loop[0]=  0;
     ind_loop[1]=  param_int[nd][IJKV];
     ind_loop[2]=  0;
     ind_loop[3]=  param_int[nd][IJKV+1];
     ind_loop[4]=  0;
     ind_loop[5]=  param_int[nd][IJKV+2];
     
      for (E_Int k = ind_loop[4]; k < ind_loop[5]; ++k)
       {
	 for (E_Int j = ind_loop[2]; j < ind_loop[3]; ++j)
  
	   {

	     i_grad = 0;
	     cflref_i = cfl[nd][j*param_int[nd][IJKV] + k*param_int[nd][IJKV]*param_int[nd][IJKV+1]];

     

	     for (E_Int i = ind_loop[0]; i < ind_loop[1]; ++i)
	       {

		 indice = i + j*param_int[nd][IJKV] + k*param_int[nd][IJKV]*param_int[nd][IJKV+1];
		 
		 cfl_loc = cfl[nd][indice];
		 if (cfl_loc > cflmax){cflmax = cfl_loc;}
		 if (cfl_loc < cflmin){cflmin = cfl_loc;}


		 /// Pour chaque couple (k,j) on cherche le nombre de cellules dans la direction i pour passer d'un niveau en temps à un autre
		 if (cfl_loc/cflref_i > 0.51 and cfl_loc/cflref_i < 1.99 ){i_grad += 1;}
           	 else
		   {
		     E_Int i_grad_min_ = min(i_grad_min,i_grad);
		     i_grad_min = i_grad_min_;
		     i_grad = 1;
		     cflref_i = cfl_loc;
		   }


				     
     

	       }	 

	   }	 

       }



     for (E_Int k = ind_loop[4]; k < ind_loop[5]; ++k)
       {
	 for (E_Int i = ind_loop[0]; i < ind_loop[1]; ++i)

	   {

	     j_grad = 0;
	     cflref_j = cfl[nd][i +  k*param_int[nd][IJKV]*param_int[nd][IJKV+1]];

	     for (E_Int j = ind_loop[2]; j < ind_loop[3]; ++j)
	       {

		 indice = i + j*param_int[nd][IJKV] + k*param_int[nd][IJKV]*param_int[nd][IJKV+1];
		 
		 cfl_loc = cfl[nd][indice];

		 
		 /// Pour chaque couple (k,i) on cherche le nombre de cellules dans la direction j pour passer d'un niveau en temps à un autre
		 if (cfl_loc/cflref_j > 0.51 and cfl_loc/cflref_j < 1.99 ){j_grad += 1;}
           	 else
		   {
		     j_grad_min = min(j_grad_min,j_grad);
		     j_grad = 1;
		     cflref_j = cfl_loc;
		   }
					     
     

	       }	 

	   }	 

       }

 
     for (E_Int j = ind_loop[2]; j < ind_loop[3]; ++j)
       {
	 for (E_Int i = ind_loop[0]; i < ind_loop[1]; ++i)

	   {

	     k_grad = 0;
	     cflref_k = cfl[nd][i +  j*param_int[nd][IJKV]];

	     for (E_Int k = ind_loop[4]; k < ind_loop[5]; ++k)
	       {

		 indice = i + j*param_int[nd][IJKV] + k*param_int[nd][IJKV]*param_int[nd][IJKV+1];
		 
		 cfl_loc = cfl[nd][indice];


		 
		 /// Pour chaque couple (j,i) on cherche le nombre de cellules dans la direction j pour passer d'un niveau en temps à un autre
		 if (cfl_loc/cflref_k > 0.50 and cfl_loc/cflref_k < 2.0 ){k_grad += 1;}
           	 else
		   {
		     k_grad_min = min(k_grad_min,k_grad);
		     k_grad = 1;
		     cflref_k = cfl_loc;
		   }


		 //if (j==0 and i == 28){ cout << cfl_loc/cfl[nd][i +  j*param_int[nd][IJKV]] << "  " << k << "  " << k_grad_min <<  endl; }
     

	       }

	     //cout << k_grad_min << "  " << j << "  " << i << endl;

	   }	 

       }



     cout << i_grad_min << "   " << j_grad_min <<"  " << k_grad_min << endl;
     //taille_bloc_i[nd] = i_grad_min;
     //taille_bloc_j[nd] = j_grad_min;
     //taille_bloc_k[nd] = k_grad_min;
     

   }

  dtmin = 1.0/cflmax;
  dtmax = 1.0/cflmin;
  exposant_max = log(dtmax / dtmin)/log(2);

  cout << "dt_min =   " << dtmin << endl ;



  

  E_Int nbbloc_i[nidom];
  E_Int nbbloc_j[nidom];
  E_Int nbbloc_k[nidom];
  E_Float cfl_max_loc;
  E_Int indice_tab; 
  E_Int maxtab = 100000000;
  E_Float dtmin_loc;


  for (E_Int nd = 0; nd <nidom ; nd++) /// Calcul nbre de blocs ds les directions i,j,k

    {

      nbbloc_i[nd] = param_int[nd][IJKV]/taille_bloc;
      if (nbbloc_i[nd]*taille_bloc < param_int[nd][IJKV] ){nbbloc_i[nd] = nbbloc_i[nd] + 1;}

      nbbloc_j[nd] = param_int[nd][IJKV+1]/taille_bloc;
      if (nbbloc_j[nd]*taille_bloc < param_int[nd][IJKV+1] ){nbbloc_j[nd] = nbbloc_j[nd] + 1;}

      nbbloc_k[nd] = param_int[nd][IJKV+2]/taille_bloc;
      if (nbbloc_k[nd]*taille_bloc < param_int[nd][IJKV+2] ){nbbloc_k[nd] = nbbloc_k[nd] + 1;}

      if (nbbloc_i[nd]*nbbloc_j[nd]*nbbloc_k[nd] > maxtab ){maxtab = nbbloc_i[nd]*nbbloc_j[nd]*nbbloc_k[nd];}


      cout << nbbloc_i[nd] << "  " << nbbloc_j[nd] << "   " << nbbloc_k[nd] << endl;


    }  





  

  E_Int tableau_level[nidom][maxtab];
  bool  tableau_bool[nidom][maxtab];
  E_Int niveau_tps;
  E_Int indice_max = 0;
  E_Int exposant;
  E_Float niv;



  for (E_Int nd = 0; nd <nidom ; nd++) /// Calcul du CFLmax et du dtmin dans chaque bloc

    {

     for (E_Int kbloc=0; kbloc < nbbloc_k[nd]; kbloc++)
        {
 	 for (E_Int jbloc=0; jbloc < nbbloc_j[nd]; jbloc++)
 	   {
 	     for (E_Int ibloc=0; ibloc < nbbloc_i[nd]; ibloc++)
 	       {
 		cfl_max_loc=0.0;

 		 for (E_Int k = kbloc*taille_bloc; k < min((kbloc+1)*taille_bloc, param_int[nd][IJKV+2]); k++)
 		   {
       		     for (E_Int j = jbloc*taille_bloc; j < min((jbloc+1)*taille_bloc, param_int[nd][IJKV+1]); j++)
 		       {	   
 			 for (E_Int i = ibloc*taille_bloc; i < min((ibloc+1)*taille_bloc, param_int[nd][IJKV]); i++)
 			 {
  			   indice = i + j*param_int[nd][IJKV] + k*param_int[nd][IJKV]*param_int[nd][IJKV+1];
 			   cfl_loc = cfl[nd][indice];
 			   if (cfl_loc > cfl_max_loc){cfl_max_loc = cfl_loc;}
			   dtmin_loc = 1.0/cfl_max_loc;
			   
 			 }
 		       }
 		   }
 
 		   indice_tab = ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd];
 		   niv = log((256*dtmin)/dtmin_loc) / log(2);
 		   if (niv < 0){niv=0.0;}
		   niveau_tps = ceil(niv);

		   //cout << dtmin_loc << "  "  << niveau_tps <<  endl;
		   
 		   tableau_level[nd][indice_tab] = niveau_tps;
 		   tableau_bool[nd][indice_tab]  = false;
		   /*
		   E_Int ideb = ibloc*taille_bloc + 1;
		   E_Int ifin = min( (ibloc+1)*taille_bloc+1 , param_int[nd][IJKV]+1);
		   E_Int jdeb = jbloc*taille_bloc + 1;
		   E_Int jfin =  min( (jbloc+1)*taille_bloc + 1 , param_int[nd][IJKV+1]+1);
		   E_Int kdeb = 1;
		   E_Int kfin = 11;

		   PyObject* l = Py_BuildValue("[i,i,i,i,i,i,i,i]",ideb,ifin,jdeb,jfin,kdeb,kfin,niveau_tps,nd);

		   PyList_Append( list_infos, l );

		   Py_DECREF( l );
		   */

 	       }

 	   }
	}
    }

  

 // Verification des rapports entre les niveaux en temps adajcents
  for (E_Int nd = 0; nd < nidom ; nd++)
    
    {

      cout << "Changement niveau sur la zone  " << nd << endl;
      
      E_Int chgt = 1;
      while (chgt != 0)

	{

	  chgt = 0;
	  
	  for (E_Int kbloc=0; kbloc < nbbloc_k[nd]; kbloc++)
	    {
	      for (E_Int jbloc=0; jbloc < nbbloc_j[nd]; jbloc++)
		{
		  for (E_Int ibloc=0; ibloc < nbbloc_i[nd]; ibloc++)
		    {
		 
		      indice = ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd];

		      
		      // zone de gauche
		      // la zone de gauche a un niveau en temps trop grand
		      // on lui donne un niveau en temps 2x plus petit que la zone courante
		      if (ibloc > 0)
			    
			{
	   
			  if ( tableau_level[nd][indice] - tableau_level[nd][indice-1] > 1)
			    
			    {
			      tableau_level[nd][indice-1] = tableau_level[nd][indice]-1;    
			      cout << "changement niveau i-1" <<  endl;
			      chgt = 1;
				
			     }
			       
			}
			 

		      // zone de droite
		      // la zone de droite a un niveau en temps trop grand
		      // on lui donne un niveau en temps 2x plus petit que la zone courante
		      if (ibloc < nbbloc_i[nd] -1) 

			{
     	  
			  if (tableau_level[nd][indice] - tableau_level[nd][indice+1] > 1) 
			    {

			      tableau_level[nd][indice + 1]=tableau_level[nd][indice] - 1; 
			      cout << "changement niveau i+1"  << endl;
			      chgt = 1;

			     }
 
			}

		      // zone bas
		      if (jbloc > 0)
			    
			{
	   
			  if ( tableau_level[nd][indice] - tableau_level[nd][indice - nbbloc_i[nd]] > 1) 
			    {
			      tableau_level[nd][indice - nbbloc_i[nd]] = tableau_level[nd][indice]-1; 
			      cout << "changement niveau j-1" << endl;
			      chgt = 1;
			    }
			       
			}
			 

		      // zone haut
		      if (jbloc < nbbloc_j[nd] -1)

			{
     	  
			  if (tableau_level[nd][indice] - tableau_level[nd][indice + nbbloc_i[nd]] > 1) 
			    {

			      tableau_level[nd][indice + nbbloc_i[nd]]=tableau_level[nd][indice] - 1; 
			      cout << "changement niveau j+1" << endl;
			      chgt = 1;
			    }
 
			}



		      // // zone arriere
		      if (kbloc > 0)
			    
		      	{
	   
		      	  if ( tableau_level[nd][indice] - tableau_level[nd][indice - nbbloc_i[nd]*nbbloc_j[nd]] > 1) 
		      	    {
		      	      tableau_level[nd][indice - nbbloc_i[nd]*nbbloc_j[nd]] = tableau_level[nd][indice]-1; 
		      	      cout << "changement niveau k-1" << endl;
		      	      chgt = 1;
		      	    }
			       
		      	}
			 

		      // zone avant
		      if (kbloc < nbbloc_k[nd] -1)

		      	{
     	  
		      	  if (tableau_level[nd][indice] - tableau_level[nd][indice + nbbloc_i[nd]*nbbloc_j[nd]] > 1) 
		      	    {

		      	      tableau_level[nd][indice + nbbloc_i[nd]*nbbloc_j[nd]]=tableau_level[nd][indice] - 1; 
		      	      cout << "changement niveau k+1" << endl;
		      	      chgt = 1;
		      	    }
 
		      	}

		      

		    }
		}
    
	    }
	}

    }


    
  
  // E_Int ind;
  // E_Int ind_depart;
  // E_Int imin;
  // E_Int imax;
  // E_Int jmin;
  // E_Int jmax;
  // E_Int kmax;
  // E_Int recul=0;


  //  for (E_Int nd = 0; nd < nidom ; nd++)

  //   {
  //      for (E_Int kbloc=0; kbloc < nbbloc_k[nd]; kbloc++)
  //       {
  //  	 for (E_Int jbloc=0; jbloc < nbbloc_j[nd]; jbloc++)
  //  	   {
  //  	     for (E_Int ibloc=0; ibloc < nbbloc_i[nd]; ibloc++) /// On parcourt tous les blocs
  //  	       {
		 
  //  		 //cout << recul << endl;
  //  		 if (recul==1){ibloc = ibloc -1;}
		 
  //  		 ind = ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd];

  //  		 recul=0;
		 
  //  		 if (tableau_bool[nd][ind]==false) /// Si le bloc courant n a pas deja ete merge on va chercher a le merger

  //  		   {

	     
  //  		     E_Int i=1;
  //  		     while (tableau_level[nd][ind] == tableau_level[nd][ind + i] and i+ibloc <= nbbloc_i[nd]-1 and tableau_bool[nd][ind + i] == false) /// Recherche du dernier bloc a merger sur la ligne i

  //  		       {
			 
  //  			 i = i+1;
		 
  //  		       }

  //  		     imax = i;
		     
  //  		     E_Int tabji[imax];
  //  		     for (E_Int i=0 ; i < imax ; i++) /// Pour chaque i on recherche le dernier bloc j a merger sur la colonne et on met j dans tabji

  //  		       {


  //  			 E_Int j=1;
  //  		     	 while (tableau_level[nd][ind] == tableau_level[nd][ind + i + j*nbbloc_i[nd]] and j+jbloc  <= nbbloc_j[nd]-1  and tableau_bool[nd][ind + i + j*nbbloc_i[nd]] == false)

  //  		     	   {
			 
  //  		     	     j = j+1; 
		 
  //  		     	   }

  //  			 tabji[i]=j;


  //  		       }

		     


  //  		     //cout << imax <<"  " << tableau_level[nd][ind]<< endl;
		     
  //  		     E_Int nbbloc_max=0;
  //  		     E_Int imax_boucle = imax;
  // 		     E_Int jmax_k=0;
  //  		     for (E_Int i=0 ; i < imax_boucle ; i++) 

  //  		       {


  //  			 E_Int jmin_k=nbbloc_j[nd];
  //  			 for (E_Int k=0; k < i+1 ; k++) // On cherche le minimum de j entre 0 et i

  //  			   {

  //  			     if(tabji[k] < jmin_k){jmin_k = tabji[k];}
  // 			     if(tabji[k] > jmax_k){jmax_k = tabji[k];}
			     
  //  			   }



  //  			 if(jmin_k*(i+1) > nbbloc_max)

  //  			   {
  //  			     recul=0;
  //  			     nbbloc_max = jmin_k*(i+1);
  //  			     imax = i+1;
  //  			     imin = 0;
  //  			     jmax = jmin_k;
  //  			     jmin=0;
			     
				 
  //  			   }

			 
			 


			 
  //  			 /*
  //  			 jmin_k=nbbloc_j[nd];
  //  			 for (E_Int k=i; k < imax_boucle ; k++) // On cherche le minimum de j entre i et imax

  //  			   {

  //  			     if(tabji[k] < jmin_k){jmin_k = tabji[k];}
			     
  //  			   }



  //  			 if(jmin_k*(imax_boucle-i) > nbbloc_max)

  //  			   {
  //  			     recul=1;
  //  			     cout << jmin_k*(imax_boucle-i) << "  " <<  nbbloc_max  << endl;
  //  			     nbbloc_max = jmin_k*(imax_boucle-i);
  //  			     imax = imax_boucle;
  //  			     imin = i;
  //  			     jmax = jmin_k;
  //  			     jmin=0;
			     
				 
  //  			   }
			 
  //  			 */
		 			
			 
  //  		       }

  // 		     E_Float rapport = float(nbbloc_max)/(float(imax)*float(jmax_k));
   		     
  // 		     printf("rapport= %4.2f \n",rapport);
		     
  //  		     cout << jmax << "  " << jmax_k << endl;
		     			 
  //  		     cout << "IMIN=  "<< imin << " IMAX= " << imax << "  JMAX= " << jmax << endl;
			 
		     

  //  	 for (E_Int j = jbloc + jmin ; j < jbloc + jmax ; j++) /// Les blocs merges sont fixes a True
  //   	     {
  //   	       for (E_Int i = ibloc + imin ; i < ibloc + imax ; i++)
  //  	        {
  //  		  tableau_bool[nd][i + j*nbbloc_i[nd]] = true ;
		    
  //  		  cout << tableau_level[nd][i + j*nbbloc_i[nd]] << endl;
  //  		 }	


  //   	    }
	     

  //  	   E_Int ideb = (ibloc+imin)*taille_bloc + 1;
  //  	   E_Int ifin = min( (ibloc+imax)*taille_bloc+1 , param_int[nd][IJKV]+1);
  //  	   E_Int jdeb = (jbloc+jmin)*taille_bloc + 1;
  //  	   E_Int jfin =  min( (jbloc+jmax)*taille_bloc + 1 , param_int[nd][IJKV+1]+1);
  //  	   E_Int kdeb = 1;
  //  	   E_Int kfin = 652;



  //  	   PyObject* l = Py_BuildValue("[i,i,i,i,i,i,i,i]",ideb,ifin,jdeb,jfin,kdeb,kfin,tableau_level[nd][ind],nd);

  //  	   PyList_Append( list_infos, l );

  //          Py_DECREF( l );
		     
			 
  //  		   }

  //  	       }

  //  	   }

  //  	}


  //   }




   //////// version avec merge en k //////////////


   
  
  E_Int ind;
  E_Int ind_depart;
  E_Int imin=0;
  E_Int imax=1;
  E_Int jmin=0;
  E_Int jmax=1;
  E_Int kmin=0;
  E_Int kmax=1;

  E_Int recul=0;


   for (E_Int nd = 0; nd < nidom ; nd++)

    {
       for (E_Int kbloc=0; kbloc < nbbloc_k[nd]; kbloc++)
        {
   	 for (E_Int jbloc=0; jbloc < nbbloc_j[nd]; jbloc++)
   	   {
   	     for (E_Int ibloc=0; ibloc < nbbloc_i[nd]; ibloc++) /// On parcourt tous les blocs
	       
   	       {
	 
		 //cout << ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd] << "   " <<  tableau_bool[nd][ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd]] << endl;

		 ind = ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd];
		 
   		 if (tableau_bool[nd][ind]==false) /// Si le bloc courant n a pas deja ete merge on va chercher a le merger


   		   {

		     
		     
   		     E_Int i=1;
   		     while (tableau_level[nd][ind] == tableau_level[nd][ind + i] and i+ibloc <= nbbloc_i[nd]-1 and tableau_bool[nd][ind + i] == false) /// Recherche du dernier bloc a merger sur la ligne i

   		       {
			 
   		 	 i = i+1;
		 
   		       }

   		     imax = i;
		     
   		     E_Int tabji[imax];
  		     for (E_Int i=0 ; i < imax ; i++) /// Pour chaque i on recherche le dernier bloc j a merger sur la colonne et on met j dans tabji

   		       {


   		 	 E_Int j=1;
   		     	 while (tableau_level[nd][ind] == tableau_level[nd][ind + i + j*nbbloc_i[nd]] and j+jbloc  <= nbbloc_j[nd]-1  and tableau_bool[nd][ind + i + j*nbbloc_i[nd]] == false)

   		     	   {
			 
   		     	     j = j+1; 
		 
   		     	   }

   		 	 tabji[i]=j;


   		       }

		     
   		     E_Int tabki[imax];
   		     for (E_Int i=0 ; i < imax ; i++)

  		       {

		 	 E_Int kmin_loc = nbbloc_k[nd];
		 	 for (E_Int j=0 ; j < tabji[i] ; j++) /// Pour chaque j on recherche le dernier bloc k

		 	   {

			 
		 	     E_Int k=1;

		 	     while (tableau_level[nd][ind] == tableau_level[nd][ind + i + j*nbbloc_i[nd] + k*nbbloc_i[nd]*nbbloc_j[nd]] and k+kbloc  <= nbbloc_k[nd]-1  and tableau_bool[nd][ind + i + j*nbbloc_i[nd] + k*nbbloc_i[nd]*nbbloc_j[nd]] == false)

		 	       {
			 
		 		 k = k+1; 
		 
		 	       }

		 	     if (k <  kmin_loc){ kmin_loc = k;} // Recherche du plus petit k parmi les j

		 	   }

		 	 tabki[i]=kmin_loc; //le plus petit k est mis dans tabki 
		       
  		       }

		     
   		     E_Int nbbloc_max=0;
   		     E_Int imax_boucle = imax;

		     
   		     E_Int jmin_loc=nbbloc_j[nd];
		     E_Int kmin_loc=nbbloc_k[nd];

   		     for (E_Int i=0 ; i < imax_boucle ; i++) // On cherche le minimum de j et le minimum de k entre 0 et i

   		       {

			 
   		 	 if(tabji[i] < jmin_loc){jmin_loc = tabji[i];}
  		 	 if(tabki[i] < kmin_loc){kmin_loc = tabki[i];}
   			   

		
   		 	 if(kmin_loc*jmin_loc*(i+1) > nbbloc_max)

   		 	   {
   		 	     recul=0;
   		 	     nbbloc_max = kmin_loc*jmin_loc*(i+1);
   		 	     imax = i+1;
   		 	     imin = 0;
   		 	     jmax = jmin_loc;
   		 	     jmin=0;
   		 	     kmax = kmin_loc;
   		 	     kmin=0;
			     
				 
   		 	   }

			 
		 
   		       }

		   
  		     //E_Float rapport = float(nbbloc_max)/(float(imax)*float(jmax_k));
   		     
  		     // printf("rapport= %4.2f \n",rapport);
		     
   		     //cout << jmax << "  " << jmax_k << endl;
		     			 
   		     //cout << "IMIN=  "<< imin << " IMAX= " << imax << "  JMIN= " << jmin << "  JMAX= " << jmax <<  "  KMIN= " << kmin << "  KMAX= " << kmax << endl;
			 
		     

		     for (E_Int k = kbloc + kmin ; k < kbloc + kmax ; k++) /// Les blocs merges sont fixes a True

		       {

		 	 for (E_Int j = jbloc + jmin ; j < jbloc + jmax ; j++) 

		 	   {

		 	     for (E_Int i = ibloc + imin ; i < ibloc + imax ; i++)

		 	       {
		   
		 		 tableau_bool[nd][i + j*nbbloc_i[nd] + k*nbbloc_i[nd]*nbbloc_j[nd]] = true ;
		    

			
		 	       }	


		 	   }


		       }


		 
		 

		     E_Int ideb = (ibloc+imin)*taille_bloc + 1;
		     E_Int ifin = min( (ibloc+imax)*taille_bloc+1 , param_int[nd][IJKV]+1);
		     E_Int jdeb = (jbloc+jmin)*taille_bloc + 1;
		     E_Int jfin =  min( (jbloc+jmax)*taille_bloc + 1 , param_int[nd][IJKV+1]+1);
		     E_Int kdeb = (kbloc+kmin)*taille_bloc + 1;
		     E_Int kfin =  min( (kbloc+kmax)*taille_bloc + 1 , param_int[nd][IJKV+2]+1);



		     PyObject* l = Py_BuildValue("[i,i,i,i,i,i,i,i]",ideb,ifin,jdeb,jfin,kdeb,kfin,tableau_level[nd][ind],nd);

		     PyList_Append( list_infos, l );

		     Py_DECREF( l );
		     
			 
		     }

   	       }

   	   }

   	}


    }





   
   
 //  E_Int ind;
 //  E_Int ind_depart;
 //  E_Int imin;
 //  E_Int imax;
 //  E_Int jmin;
 //  E_Int jmax;
 //  E_Int kmax;
 //  E_Int recul=0;

  
   
 // for (E_Int nd = 0; nd < nidom ; nd++)

 //    {
 //       for (E_Int kbloc=0; kbloc < nbbloc_k[nd]; kbloc++)
 //        {
 //  	 for (E_Int ibloc=0; ibloc < nbbloc_i[nd]; ibloc++)
 //  	   {
 //  	     for (E_Int jbloc=0; jbloc < nbbloc_j[nd]; jbloc++) /// On parcourt tous les blocs
 //  	       {
		 
 //  		 if (recul==1){jbloc = jbloc -1;}
		 
 //  		 ind = ibloc + jbloc*nbbloc_i[nd] + kbloc*nbbloc_i[nd]*nbbloc_j[nd];

 //  		 recul=0;
		 
 //  		 if (tableau_bool[nd][ind]==false) /// Si le bloc courant n a pas deja ete merge on va chercher a le merger

 //  		   {

	     
 //  		     E_Int j=1;
 //  		     while  (tableau_level[nd][ind] == tableau_level[nd][ind + j*nbbloc_i[nd]] and j+jbloc  <= nbbloc_j[nd]-1  and tableau_bool[nd][ind + j*nbbloc_i[nd]] == false) /// Recherche du dernier bloc a merger sur la ligne i

 //  		       {
			 
 //  			 j = j+1;
		 
 //  		       }

 //  		     jmax = j;
		     
 //  		     E_Int tabji[jmax];
 //  		     for (E_Int j=0 ; j < jmax ; j++) /// Pour chaque i on recherche le dernier bloc j a merger sur la colonne et on met j dans tabji

 //  		       {


 //  			 E_Int i=1;
 //  		     	 while (tableau_level[nd][ind + j*nbbloc_i[nd]] == tableau_level[nd][ind  + j*nbbloc_i[nd] + i] and i+ibloc <= nbbloc_i[nd]-1 and tableau_bool[nd][ind  + j*nbbloc_i[nd] + i] == false)

 //  		     	   {
			 
 //  		     	     i = i+1; 
		 
 //  		     	   }

 //  			 tabji[j]=i;


 //  		       }


 //  		     //cout << imax <<"  " << tableau_level[nd][ind]<< endl;
		     
 //  		     E_Int nbbloc_max=0;
 //  		     E_Int jmax_boucle = jmax;
 //  		     for (E_Int j=0 ; j < jmax_boucle ; j++) 

 //  		       {

			 
 //  			 E_Int imin_k=nbbloc_i[nd];
 //  			 for (E_Int k=0; k < j+1 ; k++) // On cherche le minimum de j entre 0 et i

 //  			   {

 //  			     if(tabji[k] < imin_k){imin_k = tabji[k];}
			     
 //  			   }



 //  			 if(imin_k*(j+1) > nbbloc_max)

 //  			   {
 //  			     recul=0;
 //  			     nbbloc_max = imin_k*(j+1);
 //  			     imax = imin_k;
 //  			     imin = 0;
 //  			     jmax = j+1;
 // 			     jmin=0;
			     
				 
 //  			   }
			 
 //  			 imin_k=nbbloc_i[nd];
 //  			 for (E_Int k=j; k < jmax_boucle ; k++) // On cherche le minimum de j entre i et imax

 //  			   {

 // 			   if(tabji[k] < imin_k){imin_k = tabji[k];}
			     
 // 			  }



 //  			 if(imin_k*(jmax_boucle-j) > nbbloc_max)

 // 			 {
 // 			   recul=1;
 // 			   cout << j << "  " << jmax_boucle   << endl;
 // 			   nbbloc_max = imin_k*(jmax_boucle-j);
 // 			   jmax = jmax_boucle;
 // 			   jmin = j;
 // 			   imax = imin_k;
 // 			   imin=0;
			     
				 
 //  			   }
			 
			 
		 			
			 
 // 			  }

		     			 
 //  		     cout << "IMIN=  "<< imin << " IMAX= " << imax << "  JMIN= " << jmin << "JMAX= " << jmax << endl;
			 
		     

 //  	 for (E_Int j = jbloc +jmin ; j < jbloc + jmax ; j++) /// Les blocs merges sont fixes a True
 //    	     {
 //    	       for (E_Int i = ibloc + imin ; i < ibloc + imax ; i++)
 //  	        {
 //  		  tableau_bool[nd][i + j*nbbloc_i[nd]] = true ;
		    
 //  		  cout << tableau_level[nd][i + j*nbbloc_i[nd]] << endl;
 //  		 }	


 //    	    }
	     

 //  	   E_Int ideb = (ibloc+imin)*taille_bloc + 1;
 //  	   E_Int ifin = min( (ibloc+imax)*taille_bloc+1 , param_int[nd][IJKV]+1);
 //  	   E_Int jdeb = (jbloc+jmin)*taille_bloc + 1;
 //  	   E_Int jfin =  min( (jbloc+jmax)*taille_bloc + 1 , param_int[nd][IJKV+1]+1);
 //  	   E_Int kdeb = 1;
 //  	   E_Int kfin = 11;



 //  	   PyObject* l = Py_BuildValue("[i,i,i,i,i,i,i,i]",ideb,ifin,jdeb,jfin,kdeb,kfin,tableau_level[nd][ind],nd);

 //  	   PyList_Append( list_infos, l );

 //           Py_DECREF( l );
		     
			 
 //  		   }

 //  	       }

 //  	   }

 //  	}


 //    }
   
  

// //cout << indice_max << endl;

    // E_Int indice_tab_d;
    // E_Int indice_tab_diag;
    // E_Int indice_tab_h;
    // E_Int max_bloc=0;
    // E_Int nb_bloc=0;
    // E_Int jmin = 10000000;
    // E_Int nb_zone = 0;
    // E_Int tableau_jmin[2000];
    // E_Int tableau_jmax[2000];
    // E_Int compteur = 0;
    // E_Int indice_tab_deb = 0;

    // //PyObject* infos;
    // //infos = new PyObject[7];
    // for (E_Int nd = 0; nd <nidom ; nd++)

    // {

    //  for (E_Int kbloc=0; kbloc < nbbloc_k; kbloc++)
    //     {
    // 	 for (E_Int jbloc=0; jbloc < nbbloc_j; jbloc++)
    // 	   {
    // 	     for (E_Int ibloc=0; ibloc < nbbloc_i; ibloc++)
    // 	       {

    // 		 ind = ibloc + jbloc*nbbloc_i + kbloc*nbbloc_i*nbbloc_j;
		 
    // 	if (tableau_bool[ind]==false)
    // 	  {

	  
    // 	    E_Int k =  ind/(nbbloc_i*nbbloc_j);
    // 	    E_Int j =  (ind - k*nbbloc_i*nbbloc_j)/nbbloc_i;
    // 	    E_Int i=  ind - j*nbbloc_i - k*nbbloc_i*nbbloc_j;
    // 	    E_Int i_=0;
    // 	    E_Int j_=0;
    // 	    E_Int nbbloc_max = 0;
    // 	    E_Int imax = 0;
    // 	    E_Int imin = 0;
    // 	    E_Int jmin = 0;
    // 	    E_Int jmax = 0;
    // 	    E_Int nbbloc = 0;
    // 	    E_Int jmin_max=100000000;

    // 	    while ( tableau_level[ind] == tableau_level[ind+i_] and i_ < nbbloc_i-1 and tableau_bool[ind+i_] == false)
    // 	      {
    // 		while ( tableau_level[ind + i_+j_*nbbloc_i] == tableau_level[ind] and j_ < nbbloc_j-1 and tableau_bool[ind+ i_+j_*nbbloc_i] == false)
    // 		  {
    // 		    j_ = j_ + 1;
    // 		  }

    // 		i_ = i_ + 1;

    // 		if (j_ < jmin_max){jmin_max = j_;}
    // 		nbbloc = i_*jmin_max;
    // 		//cout << jmin_max << endl;
    // 		if (nbbloc >= nbbloc_max)
    // 		  {
    // 		    nbbloc_max = nbbloc;
    // 		    imin = i;
    // 		    imax = i_;
    // 		    jmin = j;
    // 		    jmax = jmin_max;

    // 		  }
    // 		j_ = 0;
    // 	      }

    // 	    //cout << i << " "<< i+imax-1 <<" "<<j <<" "<< j+jmax-1 << "  " << ind <<  endl;

    // 	    for (E_Int j__ = j ;j__< j+jmax;j__++)
    // 	     {
    // 	       for (E_Int i__= i ;i__< i+imax;i__++)
    // 	    	  {
    // 		    //cout << tableau_level[i__+j__*nbbloc_i] << " "<<  i__+j__*nbbloc_i      << endl;
    // 	    	    tableau_bool[i__ + j__*nbbloc_i] = true ;
    // 	    	  }	


    // 	    }
	    
    // 	    E_Int ideb = imin*taille_tab + 1;
    // 	    E_Int ifin = min((i+imax)*taille_tab + 1,ind_loop[1]+1);
    // 	    E_Int jdeb = jmin*taille_tab + 1;
    // 	    E_Int jfin = min((j+jmax)*taille_tab + 1,ind_loop[3]+1);
    // 	    E_Int kdeb = 1;
    // 	    E_Int kfin = 201;

    // 	    nb_zone = nb_zone + 1;

    // 	    //cout << ideb << " "<< ifin <<" "<<jdeb <<" "<< jfin << "  " << tableau_level[imin + jmin*nbbloc_i] <<  endl;
    // 	    PyObject* l = Py_BuildValue("[i,i,i,i,i,i,i,i]",ideb,ifin,jdeb,jfin,kdeb,kfin,tableau_level[imin + jmin*nbbloc_i],nd);

    // 	    PyList_Append( list_infos, l );

    // 	    Py_DECREF( l );

    // 	  }




    //   }


 

    //  cout << nb_zone << endl;
    // //cout <<    tableau_jmin[2000] << " " << tableau_jmin[2000] << endl;
   

    // }

    
  return list_infos;
  
  Py_DECREF( list_infos );
  delete [] param_int;
  delete [] cfl;
  //delete [] tableau_level;
 

}






