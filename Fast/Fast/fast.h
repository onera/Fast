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
# ifndef _FAST_FAST_H_
# define _FAST_FAST_H_

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

using namespace K_FLD;

namespace K_FAST
{ 
  // Fonctions
  // =========
  PyObject* _computePT(              PyObject* self, PyObject* args);

  E_Int gsdr3( 
    E_Int**& ipt_param_int , E_Float**& ipt_param_real, 
    E_Int& nidom        , E_Int& nitrun       , E_Int& nstep    , E_Int& nstep_last, E_Int& nssiter      , E_Int& it_target , E_Int& first_it,
    E_Int& kimpli       , E_Int& lssiter_verif, E_Int& lexit_lu , E_Int& omp_mode  , E_Int& layer_mode   , E_Int& mpi       ,
    E_Int& nisdom_lu_max, E_Int& mx_nidom     , E_Int& ndimt_flt, E_Int& ndimt_grad, E_Int& threadmax_sdm, E_Int& mx_synchro, 
    E_Int& nb_pulse     ,
    E_Float& temps,
    E_Int* ipt_ijkv_sdm , 
    E_Int* ipt_ind_dm_omp       , E_Int* ipt_topology, E_Int* ipt_ind_CL, E_Int* ipt_lok, E_Int* verrou_lhs, E_Int& vartype, E_Float* timer_omp,
    E_Int* iptludic             , E_Int* iptlumax, 
    E_Int** ipt_ind_dm          , E_Int** ipt_it_lu_ssdom, E_Int** ipt_ng_pe           ,   E_Int** ipt_nfconn , E_Int** ipt_nfindex   ,
    E_Float* ipt_VectG          , E_Float* ipt_VectY     , E_Float** ipt_ssor          , E_Float** ipt_ssortmp, E_Int* ipt_ssor_size  , E_Float* ipt_drodmd,
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
    E_Float*   iptroflt      , E_Float*  iptroflt2        , E_Float*   iptgrad        , E_Float*   iptwig         , E_Float* iptstat_wig       , E_Float* iptflu,
    E_Float*   iptdrodm      , E_Float*  iptcoe           , E_Float* iptmules         , E_Float**& iptdelta        , E_Float**& iptro_res    , E_Float**& iptdrodm_trans  ,
    E_Int*&  ipt_param_int_tc, E_Float*& ipt_param_real_tc, E_Int*& ipt_linelets_int  , E_Float*& ipt_linelets_real,
    E_Int& taille_tabs       , E_Float*& stock            , E_Float*& drodmstock      , E_Float*& constk           , E_Float**  iptsrc,
    E_Float*   feq           , E_Float**& feq_transfer );
}
#endif
