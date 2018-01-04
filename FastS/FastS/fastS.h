/*    
    Copyright 2013-2018 Onera.

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

# ifndef _FASTS_FASTS_H_
# define _FASTS_FASTS_H_

# include "kcore.h"
# include "Fortran.h"

//# include "Zone.h"
using namespace K_FLD;

namespace K_FASTS
{ 
  PyObject* itt(                     PyObject* self, PyObject* args);
  PyObject* compute(                 PyObject* self, PyObject* args);
  PyObject* _computePT(              PyObject* self, PyObject* args);
  PyObject* _applyBC(                PyObject* self, PyObject* args);
  PyObject* PygetRange(              PyObject* self, PyObject* args);
  PyObject* display_ss_iteration(    PyObject* self, PyObject* args);
  PyObject* stockrecup(              PyObject* self, PyObject* args);
  PyObject* souszones_list(          PyObject* self, PyObject* args);
  PyObject* computePT_enstrophy(     PyObject* self, PyObject* args);
  PyObject* computePT_variables(     PyObject* self, PyObject* args);
  PyObject* computePT_gradient(      PyObject* self, PyObject* args);
  PyObject* computePT_velocity_ale(  PyObject* self, PyObject* args);
  PyObject* computePT_my(            PyObject* self, PyObject* args);
  PyObject* compute_effort(          PyObject* self, PyObject* args);
  PyObject* _movegrid(               PyObject* self, PyObject* args);
  PyObject* _motionlaw(              PyObject* self, PyObject* args);
  PyObject* compute_dpJ_dpW(         PyObject* self, PyObject* args);
 // PyObject* compute_RhsIterAdjoint(  PyObject* self, PyObject* args);
 // PyObject* compute_LhsIterAdjoint(  PyObject* self, PyObject* args);

  //===========
  // - Metric -
  //===========
  PyObject* metric(PyObject* self, PyObject* args);

  //===========
  // - creation grille cartesienne
  //===========
  PyObject* cartMesh(PyObject* self, PyObject* args);

  //===========
  // - init Numa (optimisation placement memoire sur DRAM)
  //===========
  PyObject* initNuma(PyObject* self, PyObject* args);

  //===========
  // - init var
  //===========
  PyObject* initVars(PyObject* self, PyObject* args);

  //==========
  // - State -
  //==========

  //=============
  // - compute -
  //=============
  // Compute t=n+1
  E_Int gsdr3(
    E_Int**& ipt_param_int , E_Float**& ipt_param_real, 
    E_Int& nidom        , E_Int& nitrun    , E_Int& nstep     , E_Int& nssiter      ,  E_Int& first_it    ,E_Int& kimpli       , E_Int& lssiter_verif,  E_Int& lexit_lu  ,  E_Int& omp_mode,
    E_Int& nisdom_lu_max, E_Int& mx_nidom  , E_Int& ndimt_flt , E_Int& threadmax_sdm, E_Int& mx_synchro   , 
    E_Int& nb_pulse     ,
    E_Float& temps,
    E_Int* ipt_ijkv_sdm       , 
    E_Int* ipt_ind_dm_omp     , E_Int* ipt_topology      ,
    E_Int* ipt_ind_sdm        , E_Int* ipt_ind_coe    , E_Int* ipt_ind_grad      ,
    E_Int* ipt_ind_CL     , E_Int* ipt_ind_CL119     , E_Int* ipt_lok       ,
    E_Int* iptludic, E_Int* iptlumax, 
    E_Int** ipt_ind_dm, E_Int** ipt_it_lu_ssdom,
    E_Float* ipt_cfl,
    E_Float**  iptx, E_Float**  ipty, E_Float** iptz,
    E_Float**  iptCellN, 
    E_Float**& iptro, E_Float**& iptro_m1, E_Float**&  iptrotmp,  E_Float**& iptro_sfd,
    E_Float**  iptmut, E_Float** iptdist,
    E_Float**  ipti, E_Float**  iptj, E_Float** iptk, E_Float** iptvol, 
    E_Float**  ipti0, E_Float**  iptj0, E_Float** iptk0,     
    E_Float**  ipti_df, E_Float**  iptj_df, E_Float** iptk_df, 
    E_Float**  iptvol_df, 
    E_Float**  iptventi, E_Float**  iptventj, E_Float** iptventk,  
    E_Float**& iptrdm,
    E_Float*   iptroflt, E_Float*  iptroflt2, E_Float*  iptwig, E_Float* iptstat_wig,
    E_Float*   iptdrodm, E_Float*  iptcoe   , E_Float* iptmules);

  //=======
  // - BC -
  //=======
  /* For applyBC */
  E_Int getDir(E_Int* ijkv,  E_Int* ind_fen);
  E_Int getRange( E_Int* ind_cgns,  E_Int* ind_fen, E_Int* param_int);


  E_Int getbcfromzone(PyObject* zone, 
                      std::vector<E_Int>& data_pos, 
                      std::vector<E_Float*>& data_bc, 
                      std::vector<E_Int*>& ind_bc ,  
                      std::vector<E_Int>& type_bc,
                      std::vector<PyArrayObject*>& hook);

  E_Int BCzone(E_Int& nd, E_Int& lrhs, E_Int& lcorner,
               E_Int* ipt_param_int, E_Float* ipt_param_real, E_Int& npass ,
               E_Int* ipt_ind_dm, E_Int* ipt_ind_dm_thread,
               E_Int* ipt_ind_CL  , E_Int* ipt_ind_CL119, E_Int* ishift_lu,
               E_Float*   iptrop  , E_Float*   ipti     , E_Float*  iptj    , E_Float* iptk   , E_Float* iptx, E_Float* ipty, E_Float* iptz,
               E_Float* iptventi  , E_Float* iptventj   , E_Float* iptventk);

}
#endif
