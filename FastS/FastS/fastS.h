/*    
    Copyright 2013-2020 Onera.

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

#ifdef _MPI
#include "CMP/include/pending_message_container.hpp"
#include "CMP/include/recv_buffer.hpp"
#include "CMP/include/send_buffer.hpp"
#include "setInterpTransfersD.h"
typedef typename CMP::PendingMsgContainer<CMP::RecvBuffer> RecvQueue;
typedef typename CMP::PendingMsgContainer<CMP::SendBuffer> SendQueue;
#endif

# include "connector.h"
# include "kcore.h"
# include "Fortran.h"

//# include "Zone.h"
using namespace K_FLD;

namespace K_FASTS
{ 
  // Structures
  // ==========
//  struct tuple_implicit_local
//  {
//     E_Int nidom_tot, lexit_lu, lssiter_verif;
//  };
  // Fonctions
  // =========
  PyObject* itt(                     PyObject* self, PyObject* args);
  PyObject* compute(                 PyObject* self, PyObject* args);
  PyObject* _computePT(              PyObject* self, PyObject* args);
  //PyObject* _computePT_(             PyObject* self, PyObject* args);
  PyObject* _stretch(                PyObject* self, PyObject* args);
  PyObject* _interpfromzone(         PyObject* self, PyObject* args);
  PyObject* _computePT_mut(          PyObject* self, PyObject* args);
  PyObject* _matvecPT(               PyObject* self, PyObject* args);
  PyObject* _applyBC(                PyObject* self, PyObject* args);
  PyObject* _applyBC_(               PyObject* self, PyObject* args);
  PyObject* PygetRange(              PyObject* self, PyObject* args);
  PyObject* display_ss_iteration(    PyObject* self, PyObject* args);
  PyObject* stockrecup(              PyObject* self, PyObject* args);
  PyObject* recup(                   PyObject* self, PyObject* args);
  PyObject* recup2(                  PyObject* self, PyObject* args);
  PyObject* recup3(                  PyObject* self, PyObject* args);
  PyObject* recup3para(              PyObject* self, PyObject* args);
  PyObject* recup3para_(             PyObject* self, PyObject* args);
  PyObject* recup3para_mpi(          PyObject* self, PyObject* args);
  PyObject* dtlocal(                 PyObject* self, PyObject* args);
  PyObject* dtlocal2(                PyObject* self, PyObject* args);
  PyObject* dtlocal2para(            PyObject* self, PyObject* args);
  PyObject* dtlocal2para_(           PyObject* self, PyObject* args);
  PyObject* dtlocal2para_mpi(        PyObject* self, PyObject* args);
  PyObject* prep_cfl(                PyObject* self, PyObject* args);
  PyObject* decoupe_maillage(        PyObject* self, PyObject* args);
  PyObject* souszones_list(          PyObject* self, PyObject* args);
  PyObject* souszones_list2(         PyObject* self, PyObject* args);
  PyObject* souszones_list3(         PyObject* self, PyObject* args);
  PyObject* distributeThreads(       PyObject* self, PyObject* args);
  PyObject* computePT_enstrophy(     PyObject* self, PyObject* args);
  PyObject* computePT_variables(     PyObject* self, PyObject* args);
  PyObject* computePT_gradient(      PyObject* self, PyObject* args);
  PyObject* computePT_velocity_ale(  PyObject* self, PyObject* args);
  PyObject* copy_velocity_ale(       PyObject* self, PyObject* args);
  PyObject* computePT_my(            PyObject* self, PyObject* args);
  PyObject* compute_effort(          PyObject* self, PyObject* args);
  PyObject* _movegrid(               PyObject* self, PyObject* args);
  //PyObject* _motionlaw(              PyObject* self, PyObject* args);
  PyObject* compute_dpJ_dpW(         PyObject* self, PyObject* args);
  PyObject* work_thread_distribution(PyObject* self, PyObject* args);
 // PyObject* compute_RhsIterAdjoint(  PyObject* self, PyObject* args);
 // PyObject* compute_LhsIterAdjoint(  PyObject* self, PyObject* args);

  //===========
  // - Metric -
  //===========
  PyObject* allocate_metric(PyObject* self, PyObject* args);
  PyObject*     init_metric(PyObject* self, PyObject* args);

  //===========
  // - SSOR -
  //===========

  PyObject* allocate_ssor(PyObject* self, PyObject* args);

  //===========
  // - Vtune -
  //===========

  PyObject* itt(PyObject* self, PyObject* args);

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

  void souszones_list_c( E_Int**& ipt_param_int, E_Float**& ipt_param_real, E_Int**& ipt_ind_dm, E_Int**& ipt_it_lu_ssdom, PyObject* work ,
                         E_Int* dtloc          , E_Int* ipt_iskip_lu      , E_Int lssiter_loc       , E_Int nidom    , 
                         E_Int nitrun          , E_Int nstep              , E_Int& nidom_tot        , E_Int& lexit_lu, E_Int& lssiter_verif);

  void distributeThreads_c( E_Int**& ipt_param_int,  E_Float**& ipt_param_real, E_Int**& ipt_ind_dm, 
                            E_Int& nidom          ,  E_Int& nssiter           , E_Int& mx_sszone   , E_Int& nstep, E_Int& nitrun, E_Int& display);

  E_Int topo_test( E_Int* topo, E_Int* nijk, E_Int& cell_tg, E_Int& lmin, E_Int& dim_i,  E_Int& dim_j, E_Int& dim_k);

  //=============
  // - compute -
  //=============
  // Compute t=n+1
  E_Int gsdr3( 
    E_Int**& ipt_param_int , E_Float**& ipt_param_real, 
    E_Int& nidom        , E_Int& nitrun       , E_Int& nstep    , E_Int& nstep_last, E_Int& nssiter , E_Int& it_target , E_Int& first_it,
    E_Int& kimpli       , E_Int& lssiter_verif, E_Int& lexit_lu , E_Int& omp_mode  , E_Int& layer_mode, E_Int& mpi,
    E_Int& nisdom_lu_max, E_Int& mx_nidom     , E_Int& ndimt_flt,
    E_Int& threadmax_sdm, E_Int& mx_synchro, 
    E_Int& nb_pulse     ,
    E_Float& temps,
    E_Int* ipt_ijkv_sdm , 
    E_Int* ipt_ind_dm_omp       , E_Int* ipt_topology, E_Int* ipt_ind_CL, E_Int* ipt_lok, E_Int* verrou_lhs, E_Int& vartype, E_Float* timer_omp,
    E_Int* iptludic             , E_Int* iptlumax, 
    E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom,
    E_Float* ipt_VectG          , E_Float* ipt_VectY     , E_Float** ipt_ssor          , E_Float** ipt_ssortmp, E_Int* ipt_ssor_size, E_Float* ipt_drodmd,
    E_Float* ipt_Hessenberg     , E_Float** iptkrylov    , E_Float** iptkrylov_transfer, E_Float* ipt_norm_kry, E_Float** ipt_gmrestmp, E_Float* ipt_givens,
    E_Float*   ipt_cfl          ,
    E_Float**  iptx             , E_Float**  ipty        , E_Float** iptz,
    E_Float**  iptCellN         , E_Float**  iptCellN_IBC, E_Int**  iptdegen,
    E_Float**& iptro, E_Float**& iptro_m1, E_Float**&  iptrotmp,  E_Float**& iptro_sfd,
    E_Float**  iptmut, E_Float*  ipt_mutd,
    E_Float**  ipti, E_Float**  iptj, E_Float** iptk, E_Float** iptvol, 
    E_Float**  ipti0, E_Float**  iptj0, E_Float** iptk0,     
    E_Float**  ipti_df, E_Float**  iptj_df, E_Float** iptk_df, 
    E_Float**  iptvol_df, 
    E_Float**  iptventi, E_Float**  iptventj, E_Float** iptventk,  
    E_Float**& iptrdm,
    E_Float*   iptroflt      , E_Float*  iptroflt2        , E_Float*   iptwig         , E_Float* iptstat_wig       ,
    E_Float*   iptdrodm      , E_Float*  iptcoe           , E_Float* iptmules         , E_Float**& iptdelta        , E_Float**& iptro_res    , E_Float**& iptdrodm_trans  ,
    E_Int*&  ipt_param_int_tc, E_Float*& ipt_param_real_tc, E_Int*& ipt_linelets_int  , E_Float*& ipt_linelets_real,
    E_Int& taille_tabs       , E_Float*& stock            , E_Float*& drodmstock      , E_Float*& constk           , E_Float**  iptsrc);


  //=============
  // - compute -
  //=============
  // Compute t=n+1
  void  matvec( 
    E_Int**& ipt_param_int , E_Float**& ipt_param_real, E_Int& no_vec_test, 
    E_Int& nidom        , E_Int& nitrun       , E_Int& nstep    , E_Int& nssiter , E_Int& it_target , E_Int& first_it,
    E_Int& kimpli       , E_Int& lssiter_verif, E_Int& lexit_lu , E_Int& omp_mode, E_Int& layer_mode, E_Int& mpi,
    E_Int& nisdom_lu_max, E_Int& mx_nidom     , E_Int& ndimt_flt,
    E_Int& threadmax_sdm, E_Int& mx_synchro, 
    E_Int& nb_pulse     ,
    E_Float& temps,
    E_Int* ipt_ijkv_sdm , 
    E_Int* ipt_ind_dm_omp       , E_Int* ipt_topology, E_Int* ipt_ind_CL, E_Int* ipt_lok, E_Int* verrou_lhs, E_Int& vartype, E_Float* timer_omp,
    E_Int* iptludic             , E_Int* iptlumax, 
    E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom,
    E_Float* ipt_VectG          , E_Float* ipt_VectY     , E_Float** ipt_ssor          , E_Float** ipt_ssortmp, E_Int* ipt_ssor_size, E_Float* ipt_drodmd,
    E_Float* ipt_Hessenberg     , E_Float** iptkrylov    , E_Float** iptkrylov_transfer, E_Float* ipt_norm_kry, E_Float** ipt_gmrestmp, E_Float* ipt_givens,
    E_Float*   ipt_cfl          ,
    E_Float**  iptx             , E_Float**  ipty        , E_Float** iptz,
    E_Float**  iptCellN         , E_Float**  iptCellN_IBC,
    E_Float**& iptro, E_Float**& iptro_m1, E_Float**&  iptrotmp,  E_Float**& iptro_sfd,
    E_Float**  iptmut, E_Float*  ipt_mutd,
    E_Float**  ipti, E_Float**  iptj, E_Float** iptk, E_Float** iptvol, 
    E_Float**  ipti0, E_Float**  iptj0, E_Float** iptk0,     
    E_Float**  ipti_df, E_Float**  iptj_df, E_Float** iptk_df, 
    E_Float**  iptvol_df, 
    E_Float**  iptventi, E_Float**  iptventj, E_Float** iptventk,  
    E_Float**& iptrdm,
    E_Float*   iptroflt, E_Float*  iptroflt2, E_Float*   iptwig, E_Float* iptstat_wig,
    E_Float*   iptdrodm, E_Float*  iptcoe   , E_Float* iptmules, E_Float**& iptdelta, E_Float**& iptro_res,
    E_Float**& iptdrodm_trans,
    E_Int*&    ipt_param_int_tc , E_Float*& ipt_param_real_tc, E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real);


  //==============================
  // - Transfer with CMP library -
  //==============================

  /* Call to transfers from FastS */

  void setInterpTransfersFastS(
  E_Float**& iptro_tmp    , E_Int& vartype             , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int     , E_Float**& param_real,
  //E_Float*& iptro_tmp    , E_Int& vartype             , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int     , E_Float**& param_real,
  E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& it_target    , E_Int& nidom            , E_Float*& ipt_timecount, E_Int& mpi           ,
  E_Int& nitcfg           , E_Int& nssiter             , E_Int& rk           , E_Int& exploc           , E_Int& numpassage );
  
  /* Transferts FastS Intra process */
  void setInterpTransfersIntra(E_Float**& ipt_ro, E_Int& vartype         , E_Int*& ipt_param_int   , E_Float*& ipt_param_real   ,
                              E_Int**& param_int, E_Float**& param_real  , E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& TypeTransfert, E_Int& nitrun, E_Int& nidom,
                              E_Int& NoTransfert, E_Float*& ipt_timecount,
                              E_Int& nitcfg     , E_Int& nssiter         , E_Int& rk, E_Int& exploc, E_Int& numpassage );

  #ifdef _MPI
  /* Transferts FastS Inter process */
  void setInterpTransfersInter(E_Float**& ipt_ro , E_Int& vartype        , E_Int*& ipt_param_int   , E_Float*& ipt_param_real   ,
                               E_Int**& param_int, E_Float**& param_real , E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& TypeTransfert, E_Int& nitrun, E_Int& nidom,
                               E_Int& NoTransfert,
                               std::pair<RecvQueue*, SendQueue*>*& pair_of_queue,
                               E_Float*& ipt_timecount,
                               E_Int& nitcfg, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& numpassage, E_Int& cpt_send_buffer );

  /* Get Transfert Inter process */
  void getTransfersInter(E_Float**& ipt_roD, E_Int**& param_int, E_Int*& param_int_tc , std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);

  /* Init Transfert Inter process */
  void init_TransferInter(std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);

  /* Delete Transfert Inter process */
  void del_TransferInter(std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);
  #endif


  //===== Distrib OMP Test
  
  E_Int topo_test( E_Int* topo, E_Int* nijk, E_Int& cells_tg, E_Int& lmin, E_Int& dim_i,  E_Int& dim_j, E_Int& dim_k);

  //=======
  // - BC -
  //=======
  /* For applyBC */

  E_Int getbcfromzone(PyObject* zone, 
                      std::vector<E_Int>& data_pos, 
                      std::vector<E_Float*>& data_bc, 
                      std::vector<E_Int*>& ind_bc ,  
                      std::vector<E_Int>& type_bc,
                      std::vector<PyArrayObject*>& hook);

  E_Int BCzone(E_Int& nd, E_Int& lrhs,  E_Int& nstep,E_Int& lcorner,
               E_Int* ipt_param_int, E_Float* ipt_param_real, E_Int& npass ,
               E_Int* ipt_ind_dm, E_Int* ipt_ind_dm_thread,
               E_Int* ipt_ind_CL  , E_Int* ipt_ind_CL119, E_Int* ipt_ind_CLgmres, E_Int* ishift_lu,
               E_Float*   iptrop  , E_Float*   ipti     , E_Float*  iptj        , E_Float* iptk   , E_Float* iptx, E_Float* ipty, E_Float* iptz,
               E_Float* iptventi  , E_Float* iptventj   , E_Float* iptventk, E_Float* iptro_gmres);

  E_Int BCzone_d(E_Int& nd, E_Int& lrhs,  E_Int& nstep,E_Int& lcorner,
               E_Int* ipt_param_int, E_Float* ipt_param_real, E_Int& npass ,
               E_Int* ipt_ind_dm, E_Int* ipt_ind_dm_thread,
               E_Int* ipt_ind_CL  , E_Int* ipt_ind_CL119, E_Int* ipt_ind_CLgmres, E_Int* ishift_lu,
               E_Float*   iptrop  , E_Float*   ipti     , E_Float*  iptj        , E_Float* iptk   , E_Float* iptx, E_Float* ipty, E_Float* iptz,
               E_Float* iptventi  , E_Float* iptventj   , E_Float* iptventk, E_Float* iptro_gmres);

  void __applyBC(E_Float**& iptro,E_Int& nidom, E_Int**& ipt_param_int, E_Float**& ipt_param_real,
                 E_Float** iptx,  E_Float** ipty, E_Float** iptz,
                 E_Float** ipti,  E_Float** iptj, E_Float** iptk, E_Float** iptvol,
                 E_Float** iptventi        , E_Float**  iptventj         , E_Float** iptventk);


  //====================
  // - explicite local -
  //====================

  //void dtlocal2para_c(E_Float**& iptro, E_Float**& iptro_p1, E_Int*& param_int,
  //	    E_Float*& param_real, E_Int**& param_intt, E_Float**& param_realt,
  //		    E_Float*& iptdrodm,  E_Float*& iptdrodmcoe, E_Float*& iptstk, E_Float*& iptdrodmstk,
  //		    E_Float*& iptcstk, E_Int& nstep, E_Int& omp_mode, E_Int& taille_tabs, E_Int& nidom, E_Int& mx_sszone, E_Int*& ipt_ssor_size);

void dtlocal2para_c(E_Float**& iptro, E_Float**& iptro_p1, E_Int*& param_int, E_Float*& param_real,  E_Int**& param_intt, E_Float**& param_realt,
                      E_Float*& iptdrodm,  E_Float*& iptdrodmcoe, E_Float*& iptstk, E_Float*& iptdrodmstk, E_Float*& iptcstk,  E_Int& nstep, 
  		      E_Int& omp_mode, E_Int& taille_tabs, E_Int& nidom);




void BC_local(E_Float**& iptro, E_Float**& iptro_p1, E_Int*& param_int,
		    E_Float*& param_real, E_Int**& param_intt, E_Float**& param_realt,
		    E_Float*& iptdrodm,  E_Float*& iptdrodmcoe, E_Int* ipt_ind_CL, E_Float** ipti, E_Float** iptj, E_Float** iptk,
		    E_Float** iptx, E_Float** ipty, E_Float** iptz, E_Float** iptventi, E_Float** iptventj, E_Float** iptventk, 
		    E_Float*& iptstk, E_Float*& iptdrodmstk,
		    E_Float*& iptcstk, E_Int& nstep, E_Int& omp_mode, E_Int& taille_tabs, E_Int& nidom);

void recup3para_c(E_Float**& iptro, E_Int*& param_int,
		    E_Float*& param_real, E_Int**& param_intt, 
		    E_Float*& iptstk, E_Int& nstep, E_Int& omp_mode, E_Int& taille_tabs, E_Int& nidom);



}
#endif
