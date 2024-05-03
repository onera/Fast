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
/*# include "fastS.h"*/
# include "fast.h"
//# include "../FastASLBM/FastASLBM/fastLBM.h"
# include "param_solver.h"
# include <string.h>
# include <omp.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
PyObject* K_FAST::filtering_(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* work; PyObject* isFilter;
  E_Int loc, nstep, nitrun, omp_mode;

#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOlll" , &zones , &work, &isFilter, &nstep, &nitrun, &omp_mode)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "OOOiii" , &zones , &work, &isFilter, &nstep, &nitrun, &omp_mode)) return NULL;
#endif


  //Nombre de zones
  E_Int nidom  = PyList_Size(zones);
  cout << nidom << endl;

  E_Int**   param_intt;

  E_Float** param_realt;  E_Float** iptro_p1;   E_Float** iptro; E_Float** iptFilterN;

  param_intt  = new E_Int*[nidom];

  param_realt = new E_Float*[nidom*4];
  iptro       = param_realt + nidom;
  iptro_p1    = iptro       + nidom;
  iptFilterN  = iptro_p1    + nidom;

  vector<PyArrayObject*> hook;

  E_Int neq_max       = 5;
  E_Int flag_filter   = 0;
  E_Int ndim_filter   = 0;
  E_Int ighost_filter = 0;

  // On parcourt les zones et on recupere les tableaux/vecteurs
  for (E_Int nd = 0; nd < nidom; nd++)
     {
       // check zone
       PyObject* zone = PyList_GetItem(zones, nd);

       /* Get numerics from zones */
       PyObject* numerics   = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*       o    = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"  );
       param_intt[nd]       = K_PYTREE::getValueAI(o, hook);
	                     o    = K_PYTREE::getNodeFromName1(numerics, "Parameter_real" );
       param_realt[nd]      = K_PYTREE::getValueAF(o, hook);

       /*Pointeur var primitives */
       PyObject* sol_center = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
       PyObject*         t  = K_PYTREE::getNodeFromName1(sol_center, "Density");
                  iptro[nd] = K_PYTREE::getValueAF(t, hook);

       /*Pointeur idx filtrage */
       t            = K_PYTREE::getNodeFromName1(sol_center, "FilterN");
       if (t == NULL)
       {
         iptFilterN[nd] = NULL;
         flag_filter = flag_filter + 0;
       }
       else
       {
         iptFilterN[nd] = K_PYTREE::getValueAF(t, hook);
         flag_filter = flag_filter + 1;
         //Si zone a filter, on ajoute la taille de la zone aux tableaux de travail
         cout << (param_intt[nd][ NIJK ]*param_intt[nd][ NIJK+1 ]*param_intt[nd][ NIJK+2 ]) << endl;
         ndim_filter = ndim_filter + (param_intt[nd][ NIJK ]*param_intt[nd][ NIJK+1 ]*param_intt[nd][ NIJK+2 ]);
       }


     } //fin boucle zones


  //Reservation tableau travail temporaire
  PyObject* drodm2Array = PyDict_GetItemString(work,"hors_eq"); FldArrayF* drodm2;
  K_NUMPY::getFromNumpyArray(drodm2Array, drodm2, true); E_Float* iptroflt = drodm2->begin();

  E_Int shift_dom   = 0;
  E_Int shift_roflt = 0;
  flag_filter = 1;

  //cout << "Flag filter = " << flag_filter << endl;
  if (flag_filter!=0)
  {
    E_Int filtered = 0;

    //----------------------------------------------------------------------------
    //                     BOUCLE SUR LES DOMAINES PRESENTS
    //----------------------------------------------------------------------------
    for  (E_Int nd=0; nd < nidom; nd++)
        {
          cout << "Domaine = " << nd << endl;

          PyObject* my_val = PyList_GetItem(isFilter, nd);
          E_Int to_filter = (E_Int) PyInt_AsLong(my_val);
          //cout << my_val << endl;

          E_Int taille_zone = (param_intt[nd][ NIJK ]*param_intt[nd][ NIJK+1 ]*param_intt[nd][ NIJK+2 ]);

          if (to_filter==1)
          {

            E_Int nx = param_intt[nd][ IJKV   ];
          	E_Int ny = param_intt[nd][ IJKV+1 ];
          	E_Int nz = param_intt[nd][ IJKV+2 ];
          	E_Int ng = param_intt[nd][ NIJK+3 ];

            cout << nx << " " << ny << " " << nz << " " << ng << endl;

            move_to_temp_( nd, nidom, param_intt[nd], iptro[nd], iptroflt+shift_roflt); //, ng, nx, ny, nz);
            filtrage_5_(   nd, nidom, param_intt[nd], iptro[nd], iptroflt+shift_roflt); //, ng, nx, ny, nz);

            //filtered = filtered + 1;
            //shift_roflt = shift_roflt + taille_zone*5;

          }

          shift_roflt = shift_roflt + param_intt[nd][ NDIMDX ]*param_intt[nd][ NEQ ];
          //Tout ce qui suit s'applique uniquement aux zones a filter
          //cout << "On filtre ? " << flag_filter[nd] << endl;
          // if (flag_filter[nd]==1)
          // {
          //
          //   E_Int taille_zone = param_intt[nd][ NIJK ]*param_intt[nd][ NIJK+1 ]*param_intt[nd][ NIJK+2 ];
          //
          // }

          //Copie des valeurs dans les tableaux de travail
          //cp_values_filter_()

          //filtrage_(shift_roflt)
          /*indice_boucle_lu_(nd, ithread_loc, Nbre_thread_actif_loc, param_intt[nd][ ITYPCP ],
                            donorPts_, topology, ipt_ind_dm_omp_thread);

    		  copy_valuespara_( param_intt[nd], donorPts, donorPts_, tab1, tab2, ind, sol, taillefenetre);

          shift_dom = shift_dom + nb_rac_zone*taille_rac_max;
          */
        }// boucle raccords
      }

 delete [] param_intt; delete [] param_realt;

 RELEASESHAREDN( drodm2Array    , drodm2);
 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

 Py_INCREF(Py_None);
 return Py_None;
}
