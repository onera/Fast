/*    
    Copyright 2013-2022 Onera.

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
# include "fastS.h"
# include "param_solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
using namespace K_FLD;

//=============================================================================
/* 
   Applique les conditions aux limites d'un type donne sur la zone
   *in place* 
*/
//=============================================================================
PyObject* K_FASTS::compute_dpJ_dpW(PyObject* self, PyObject* args)
{
  PyObject *zones; PyObject *zones_eff; PyObject *metrics; PyObject* work; E_Float cosAoA; E_Float sinAoA; E_Float surfinv;
  if (!PyArg_ParseTuple(args, "OOOOddd", &zones, &zones_eff, &metrics, &work, &cosAoA, &sinAoA, &surfinv)) return NULL;



  E_Int**   ipt_param_int;
  
  E_Float** ipt_param_real  ;
  E_Float** iptx;       E_Float** ipty;     E_Float** iptz;    E_Float** iptro; E_Float** iptmut;
  E_Float** ipti;       E_Float** iptj;     E_Float** iptk;    E_Float** iptvol;
  E_Float** ipti0;      E_Float** iptj0;    E_Float** iptk0;
  E_Float** ipti_df;    E_Float** iptj_df;  E_Float** iptk_df ; E_Float** iptvol_df;
  E_Float** iptdist;    E_Float** iptventi; E_Float** iptventj; E_Float** iptventk;
  E_Float** ipt_dpCdp_dpW; E_Float** ipt_dpClp_dpW;

  E_Int nidom    = PyList_Size(zones);

  ipt_param_int     = new E_Int*[nidom];

  //iptx              = new E_Float*[nidom*20];
  // redimensionner en fonction du nombre de tableau a stocker
  iptx              = new E_Float*[nidom*22];
  ipty              = iptx            + nidom;
  iptz              = ipty            + nidom;
  iptro             = iptz            + nidom;
  iptmut            = iptro           + nidom;
  ipti              = iptmut          + nidom;
  iptj              = ipti            + nidom;
  iptk              = iptj            + nidom;
  iptvol            = iptk            + nidom;
  ipti0             = iptvol          + nidom;
  iptj0             = ipti0           + nidom;
  iptk0             = iptj0           + nidom;
  ipti_df           = iptk0           + nidom;
  iptj_df           = ipti_df         + nidom;
  iptk_df           = iptj_df         + nidom;
  iptvol_df         = iptk_df         + nidom;
  iptventi          = iptvol_df       + nidom;
  iptventj          = iptventi        + nidom;
  iptventk          = iptventj        + nidom;
  ipt_param_real    = iptventk        + nidom;
  ipt_dpCdp_dpW     = ipt_param_real  + nidom;
  ipt_dpClp_dpW     = ipt_dpCdp_dpW   + nidom;


  /*------*/
  /* Zone */
  /*------*/
  vector<PyArrayObject*> hook;


  /// Tableau pour stockage senseur oscillation
  PyObject* wigArray = PyDict_GetItemString(work,"wiggle"); FldArrayF* wig; FldArrayF* eff;
  K_NUMPY::getFromNumpyArray(wigArray, wig, true); E_Float* iptwig = wig->begin();
//  K_NUMPY::getFromNumpyArray(effarray, eff, true); E_Float* ipteff = eff->begin();

  
  ///
  ///
  /// Recuperation Pointeur arbre NS
  ///
  ///
  for (E_Int nd = 0; nd < nidom; nd++)
  { 
     // check zone
     PyObject* zone   = PyList_GetItem(zones  , nd); 
     PyObject* metric = PyList_GetItem(metrics, nd);


    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);
                       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(o, hook);

    //
    //
    //Pointeur maillage
    //
    //
    if(ipt_param_int[nd][ LALE ]== 0){ GET_XYZ( "GridCoordinates"     , zone, iptx[nd], ipty[nd], iptz[nd])}
    else                             { GET_XYZ( "GridCoordinates#Init", zone, iptx[nd], ipty[nd], iptz[nd])}

    //
    //Pointeur var primitive + visco + distance paroi + cellN
    //
    //
    PyObject* sol_center; PyObject* t; PyObject* sol_vertex;
    sol_center    = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
    sol_vertex    = K_PYTREE::getNodeFromName1(zone      , "FlowSolution");
    t             = K_PYTREE::getNodeFromName1(sol_center, "Density");
    iptro[nd]     = K_PYTREE::getValueAF(t, hook);
    t             = K_PYTREE::getNodeFromName1(sol_center, "dpCDp_dpDensity");
    ipt_dpCdp_dpW[nd]= K_PYTREE::getValueAF(t, hook);
    t             = K_PYTREE::getNodeFromName1(sol_center, "dpCLp_dpDensity");
    ipt_dpClp_dpW[nd]= K_PYTREE::getValueAF(t, hook);


    if(ipt_param_int[nd][ IFLOW ] > 1)
      { t  = K_PYTREE::getNodeFromName1(sol_center, "ViscosityEddy");
        if (t == NULL) { PyErr_SetString(PyExc_ValueError, "ViscosityEddy is missing for NS computations."); return NULL; }
        else           { iptmut[nd]   = K_PYTREE::getValueAF(t, hook);}
      }
    else {iptmut[nd] = iptro[nd];}

    //
    //
    //Pointeur metric
    //
    //
    GET_TI(   METRIC_TI  , metric, ipt_param_int[nd], ipti [nd]  , iptj [nd]  , iptk [nd]  , iptvol[nd]   )
    GET_TI(   METRIC_TI0 , metric, ipt_param_int[nd], ipti0[nd]  , iptj0[nd]  , iptk0[nd]  , E_Float* tmp )
    GET_TI(   METRIC_TIDF, metric, ipt_param_int[nd], ipti_df[nd], iptj_df[nd], iptk_df[nd], iptvol_df[nd] )

    GET_VENT( METRIC_VENT, metric, ipt_param_int[nd], iptventi[nd], iptventj[nd], iptventk[nd] )

  } // boucle zone arbre NS

  ///
  ///
  /// Recuperation Pointeur arbre effort
  ///
  ///
  E_Int**  ipt_param_int_eff;

  nidom             = PyList_Size(zones_eff);
  ipt_param_int_eff = new E_Int*[nidom];

  for (E_Int nd = 0; nd < nidom; nd++)
  { 

     PyObject* zone   = PyList_GetItem(zones_eff , nd); 

     /* Get numerics from zone */
     PyObject* numerics    = K_PYTREE::getNodeFromName1(zone, ".Solver#ownData");
     PyObject*        t    = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
     ipt_param_int_eff[nd] = K_PYTREE::getValueAI(t, hook);

   } // boucle zone arbre effort


  FldArrayF  xyz_ref(3); E_Float* ipt_xyz_ref = xyz_ref.begin();
  xyz_ref[0]=0.; xyz_ref[1]=0.;xyz_ref[2]=0.;


  E_Int  Nbre_thread_max = 1;
#ifdef _OPENMP
  Nbre_thread_max = omp_get_max_threads();
#endif
//  FldArrayF          effort(7*Nbre_thread_max); 
  FldArrayI thread_topology(3*Nbre_thread_max); 
  FldArrayI   ind_dm_thread(6*Nbre_thread_max);  

  E_Int lmin=4;
  //
  //
  //loop sur les fenetres pour calcul flux
  //
  //

#pragma omp parallel default(shared)
  {

#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

       for (E_Int nd = 0; nd < nidom; nd++)
       { 

       E_Int nd_ns  = ipt_param_int_eff[nd][EFF_NONZ];

       //printf("no zone ns %d %d \n", nd, nd_ns);
       E_Int shift_wig = 0;
       for (E_Int i = 0; i < nd_ns; i++) shift_wig  = shift_wig  + ipt_param_int[i][ NDIMDX ]*3;

       E_Int* ipt_thread_topology = thread_topology.begin() + (ithread-1)*3;
       E_Int* ind_loop            = ind_dm_thread.begin()   + (ithread-1)*6;
      // E_Float* effort_omp        = effort.begin()          + (ithread-1)*7;      A remplacer effort_omp par dclp et dcdp ?


      // if(nd==0) for (E_Int i = 0; i < 7; i++) effort_omp[i]=0;

        //calcul du sous domaine a traiter par le thread 

        
        indice_boucle_lu_(nd, ithread, Nbre_thread_actif, lmin,
                           ipt_param_int_eff[nd]+EFF_LOOP,
                           ipt_thread_topology, ind_loop);

        dpj_dpw_calc_meth_( nd, ithread, ipt_param_int[nd_ns], ipt_param_real[nd_ns], ipt_param_int_eff[nd],
                      ind_loop            , ipt_xyz_ref         ,
                      iptro[nd_ns]        , iptwig +shift_wig   , iptmut[nd_ns]   ,   
                      iptx[nd_ns]         , ipty[nd_ns]         , iptz[nd_ns]     ,
                      ipt_dpCdp_dpW[nd_ns], ipt_dpClp_dpW[nd_ns],
                      ipti[nd_ns]         , iptj[nd_ns]         , iptk[nd_ns]     , iptvol[nd_ns] ,
                      iptventi[nd_ns]     , iptventj[nd_ns]     , iptventk[nd_ns] , cosAoA, sinAoA, surfinv);


          // REMARQUE : COMPLEXITE AVEC LA ROUTINE TAPENADE POUR LES CALCULS PARALLELES
          //            ==> AJOUTER UN ARGUMENT A RHS_ADJOINT POUR QUE LANCER LA ROUTINE GENEREE PAR TAPENADE POUR EFFECTUER DES CALCULS.
          //            ==> PRENDRE EN COMPTE TOUS LES PARAMETRES NECESSAIRES AU CALCUL POUR LA ROUTINE TAPENADE.


       } // fin loop zone
  } // Fin zone // omp


        



  delete [] iptx; delete [] ipt_param_int; delete [] ipt_param_int_eff;

  RELEASESHAREDN( wigArray  , wig  );


  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
