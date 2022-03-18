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

# ifndef _FASTC_FASTC_H_
# define _FASTC_FASTC_H_


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

namespace K_FASTC
{ 

  // Check a value in Numerics Dictionary
  E_Int checkNumericsValue(PyObject* numerics, const char* value,
                           E_Int& retInt, E_Float& retFloat, char*& retChar);


  PyObject* _motionlaw( PyObject* self, PyObject* args);
  PyObject* PygetRange( PyObject* self, PyObject* args);


  //=======
  // - BC -
  //=======
  /* For BCcompact */

  E_Int getRange( E_Int* ind_cgns,  E_Int* ind_fen, E_Int* param_int);
  E_Int getDir(E_Int* ijkv,  E_Int* ind_fen);

  //==============================
  // - Transfer with CMP library -
  //==============================

  /* Call to transfers from FastS */

  void setInterpTransfersFast(
  E_Float**& iptro_tmp    , E_Int& vartype             , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int     , E_Float**& param_real, E_Int*& ipt_omp,
  E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& it_target    , E_Int& nidom            , E_Float*& ipt_timecount, E_Int& mpi           ,
  E_Int& nitcfg           , E_Int& nssiter             , E_Int& rk           , E_Int& exploc           , E_Int& numpassage );
  
  /* Transferts FastS Intra process */
  void setInterpTransfersIntra(E_Float**& ipt_ro, E_Int& vartype         , E_Int*& ipt_param_int   , E_Float*& ipt_param_real, 
                              E_Int**& param_int, E_Float**& param_real  , E_Int*& ipt_omp         , E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& TypeTransfert, E_Int& nitrun, E_Int& nidom,
                              E_Int& NoTransfert, E_Float*& ipt_timecount,
                              E_Int& nitcfg     , E_Int& nssiter         , E_Int& rk, E_Int& exploc, E_Int& numpassage );

  /*---
  //LBM
  ---*/

  /* Transfers FastLBM */
  void setInterpTransfersFastLBM(E_Float**& iptro_tmp       , E_Int& vartype             ,
				 E_Int*& param_int_tc       , E_Float*& param_real_tc    ,
				 E_Int**& param_int         , E_Float**& param_real      ,
				 E_Int*& ipt_linelets_int   , E_Float*& ipt_linelets_real,
				 E_Int& it_target           , E_Int& nidom               , 
				 E_Float*& ipt_timecount    , E_Int& mpi                 ,
				 E_Int& nitcfg              , E_Int& nssiter             ,
				 E_Int& rk                  , E_Int& exploc              ,
				 E_Int& numpassage          ,
				 E_Float**& ipt_macros_local, E_Float**& ipt_Qneq_local  );
  
  /* Intra process */
  void setInterpTransfersIntraLBM(E_Float**& ipt_ro          , E_Int& vartype             ,
				  E_Int*& ipt_param_int      , E_Float*& ipt_param_real   ,
				  E_Int**& param_int         , E_Float**& param_real      ,
				  E_Int*& ipt_linelets_int   , E_Float*& ipt_linelets_real,
				  E_Int& TypeTransfert       , E_Int& nitrun              ,
				  E_Int& nidom               , E_Int& NoTransfert         ,
				  E_Float*& ipt_timecount    , E_Int& nitcfg              ,
				  E_Int& nssiter             , E_Int& rk                  ,
				  E_Int& exploc              , E_Int& numpassage          ,
				  E_Float**& ipt_macros_local, E_Float**& ipt_Qneq_local  );

  #ifdef _MPI
  /* Transferts FastS Inter process */
  void setInterpTransfersInter(E_Float**& ipt_ro , E_Int& vartype        , E_Int*& ipt_param_int   , E_Float*& ipt_param_real ,
                               E_Int**& param_int, E_Float**& param_real , E_Int*& ipt_omp, E_Int*& ipt_linelets_int, E_Float*& ipt_linelets_real, E_Int& TypeTransfert, E_Int& nitrun, E_Int& nidom,
                               E_Int& NoTransfert,
                               std::pair<RecvQueue*, SendQueue*>*& pair_of_queue,
                               E_Float*& ipt_timecount,
                               E_Int& nitcfg, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& numpassage , E_Int& nb_send_buffer);

  /* Get Transfert Inter process */
  void getTransfersInter(E_Float**& ipt_roD, E_Int**& param_int, E_Int*& param_int_tc , std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);

  /* Init Transfert Inter process */
  void init_TransferInter(std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);

  /* Delete Transfert Inter process */
  void del_TransferInter(std::pair<RecvQueue*, SendQueue*>*& pair_of_queue);
  #endif

}
#endif
