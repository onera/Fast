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


#include "FastC/fastc.h"
#include "FastC/param_solver.h"

#ifdef _MPI
#include <mpi.h>
#include "CMP/include/pending_message_container.h"
#include "CMP/include/recv_buffer.hpp"
#include "CMP/include/send_buffer.hpp"
#include "TRANSFERT/setInterpTransfersD.h"
#endif


#include <utility>

using namespace std;
using namespace K_FLD;

//#define TimeShow

#ifdef TimeShow

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

E_Int timeShowFastS = 1;
#endif

static E_Int timelevel_tc=0;

#ifdef _MPI

#if defined(PROFILE_TRANSFERT)
static double send_buffer_time = 0.;
static double recv_buffer_time = 0.;
static double isend_time = 0.;
#endif

static std::pair<RecvQueue*, SendQueue*>* pair_of_queue_pass1 = NULL;
static std::pair<RecvQueue*, SendQueue*>* pair_of_queue_pass2 = NULL;
static std::pair<RecvQueue*, SendQueue*>* pair_of_queue_pass3 = NULL;
static std::pair<RecvQueue*, SendQueue*>* pair_of_queue_pass4 = NULL;

void K_FASTC::init_TransferInter(
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue_loc) {
  int szCom;
  MPI_Comm_size(MPI_COMM_WORLD, &szCom);
  RecvQueue* pt_RecvQueue = new RecvQueue(szCom);
  SendQueue* pt_SendQueue = new SendQueue(szCom);
  pair_of_queue_loc = new std::pair<RecvQueue*, SendQueue*>(pt_RecvQueue, pt_SendQueue);
//      std::make_pair(pt_RecvQueue, pt_SendQueue));
}

void K_FASTC::del_TransferInter(
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue_loc) {
  delete pair_of_queue_loc->first;
  delete pair_of_queue_loc->second;
  delete pair_of_queue_loc;
  pair_of_queue_loc = NULL;
}
#endif

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
void K_FASTC::setInterpTransfersFast(
  E_Float**& iptro_tmp, E_Int& vartype         , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int  , E_Float**& param_real  , E_Int*& ipt_omp,
  E_Int*& linelets_int, E_Float*& linelets_real, E_Int& it_target    , E_Int& nidom            , E_Float*& ipt_timecount, E_Int& mpi,
  E_Int& nstep        , E_Int& nitmax          , E_Int& rk           , E_Int& exploc           , E_Int& numpassage, E_Int& Nopass)

{
  int rank = 0;
  E_Int dest = 0;

#ifdef _MPI
  if(mpi) { MPI_Comm_rank (MPI_COMM_WORLD, &rank); }
#endif

#ifdef TimeShow
   ofstream outputfile;
   std::ostringstream tmp;
   tmp << "Output" << std::setw(4) << std::setfill('0') << std::to_string(rank);
   std::string filename = tmp.str();
   outputfile.open(filename, ios::app);
#ifdef _OPENMP
   E_Float time_in = omp_get_wtime();
#endif
#endif

 //Swap (call to setInterpTransfer)
  //if ( (param_int_tc != NULL) && (param_real_tc != NULL))
  if  (param_int_tc != NULL) 
  {
    E_Int TypeTransfert ;    
    E_Int Nbp2p_send = param_int_tc[1];
    E_Int sizecomID  = param_int_tc[2];
    E_Int shift_graph = sizecomID + 2;

    timelevel_tc=0;
    if(Nbp2p_send !=0) // evite depassememnt tableau si pas d'envoi dans param_int, seulememnt reception
     {
        E_Int ech        = param_int_tc[1 + shift_graph];
        timelevel_tc     = param_int_tc[ech + 3];
     }
    //printf("VERIF: szCom et Nbp2p_send  %d %d init= %d  %p, nstep = %d , mpi= %d , timelevel_tc %d \n", sizecomID, Nbp2p_send , param_int_tc[0], param_real_tc,  nstep, mpi, timelevel_tc);fflush(0);

    #ifdef _MPI
    //std::pair<RecvQueue*, SendQueue*>* pair_of_queue;
    RecvQueue* pt_rcv_queue     = NULL;

    E_Int etiquette=404+Nopass;
    /*
    std::pair<RecvQueue*, SendQueue*>* pair_of_queue_loc;
    if (mpi)
      { if     (Nopass==0) { pair_of_queue_loc = pair_of_queue_pass1;}
        else if(Nopass==1) { pair_of_queue_loc = pair_of_queue_pass2;}
        else if(Nopass==2) { pair_of_queue_loc = pair_of_queue_pass3;}
        else if(Nopass==3) { pair_of_queue_loc = pair_of_queue_pass4;}
        else   { printf("Error pair of queue too small: transfert failed \n");exit(0);}
      }
    */

    //premier passage dans transfert couche C depuis  mise a plat de tc
    if (param_int_tc[0]==0 and mpi)
    {  
        if     (Nopass==0 && pair_of_queue_pass1 != NULL ){ K_FASTC::del_TransferInter(pair_of_queue_pass1);}
        else if(Nopass==1 && pair_of_queue_pass2 != NULL ){ K_FASTC::del_TransferInter(pair_of_queue_pass2);}
        else if(Nopass==2 && pair_of_queue_pass3 != NULL ){ K_FASTC::del_TransferInter(pair_of_queue_pass3);}
        else if(Nopass==3 && pair_of_queue_pass4 != NULL ){ K_FASTC::del_TransferInter(pair_of_queue_pass4);}

        if (Nopass >=4)  { printf("Error pair of queue too small: del transfert failed \n");exit(0);}

        //if (pair_of_queue_loc  != NULL ) { K_FASTC::del_TransferInter(pair_of_queue_loc);}

        //printf("ALLO %d init= %d , nstep = %d , mpi= %d , NoPass: %d \n", sizecomID, param_int_tc[0], nstep, mpi, Nopass); fflush(0);
        //flag transfer initilisé. Remise a zero dans miseAplat.
        param_int_tc[0]=1;
    }

    E_Int nbcomID_S; E_Int nbcomID_U;  E_Int nbcomIBC_S; E_Int pt_debID_S; E_Int pt_debID_U; E_Int pt_debIBC_S;

    // info Comm ID instationnaire et dtloc G Jeanmass en reception
    E_Int iter = 1;
    if (exploc==1 and mpi) { iter = nstep; nbcomID_U =0;}
    else if(timelevel_tc !=0)
       {  
          pt_debID_U = param_int_tc[2+ iter + 1 + it_target] + 2 +1;
          nbcomID_U  = param_int_tc[ pt_debID_U ];
       }
    else { nbcomID_U =0;}

    // info Comm ID stationnaire en reception
    pt_debID_S = param_int_tc[2+iter] + 2 +1;
    nbcomID_S  = param_int_tc[ pt_debID_S ];


    E_Int send_mpi=0;  // flag pour connaitre existence envoi mpi
    //printf("Nb source  %d %d , pt: %d , nstep: %d , Nopass %d \n", nbcomID_S, nbcomID_U, pt_debID_S, nstep, Nopass ); fflush(0);
    if (mpi and (nbcomID_S != 0 or nbcomID_U != 0 or Nbp2p_send !=0))
    {
      if     (Nopass==0 && pair_of_queue_pass1 == NULL ){ K_FASTC::init_TransferInter(pair_of_queue_pass1); }
      else if(Nopass==1 && pair_of_queue_pass2 == NULL ){ K_FASTC::init_TransferInter(pair_of_queue_pass2); }
      else if(Nopass==2 && pair_of_queue_pass3 == NULL ){ K_FASTC::init_TransferInter(pair_of_queue_pass3); }
      else if(Nopass==3 && pair_of_queue_pass4 == NULL ){ K_FASTC::init_TransferInter(pair_of_queue_pass4); }
      if (Nopass >=4)  { printf("Error pair of queue too small: init  transfert failed \n");exit(0);}

      //if (pair_of_queue_loc == NULL ) { K_FASTC::init_TransferInter(pair_of_queue_loc );}

      #ifdef TimeShow
      #ifdef _OPENMP
       time_in = omp_get_wtime();
      #endif
      #endif
      //
      //  MPI_Barrier(MPI_COMM_WORLD);
      //
      //Debut Transfert
      //
      //
      if     (Nopass==0 ){ pt_rcv_queue = pair_of_queue_pass1->first;}
      else if(Nopass==1 ){ pt_rcv_queue = pair_of_queue_pass2->first;}
      else if(Nopass==2 ){ pt_rcv_queue = pair_of_queue_pass3->first;}
      else if(Nopass==3 ){ pt_rcv_queue = pair_of_queue_pass4->first;}
      //pt_rcv_queue = pair_of_queue_loc->first;

      if (pt_rcv_queue->size() == 0 )
        {
          for (E_Int ircv = 1; ircv < nbcomID_S +1; ++ircv)
           {
            E_Int source = param_int_tc[ pt_debID_S + ircv];
            pt_rcv_queue->emplace_back( source , etiquette);
            CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
            recv_buffer.irecv();
            //bool flag = recv_buffer.test();
            //printf("reception ID Steady source  %d %d \n", source, nstep ); fflush(0);
           }
          for (E_Int ircv = 1; ircv < nbcomID_U +1; ++ircv)
           {
            E_Int source = param_int_tc[ pt_debID_U + ircv];
            pt_rcv_queue->emplace_back( source , etiquette);
            CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
            recv_buffer.irecv();
            //bool flag = recv_buffer.test();
            //printf("reception ID Unsteady source  %d %d \n", source, nstep ); fflush(0);
           }
        }
      else
        {   pt_rcv_queue->resize(nbcomID_S);
            assert(pt_rcv_queue->size() == nbcomID_S );
            for ( auto iterBuf = pt_rcv_queue->begin(); iterBuf != pt_rcv_queue->end(); ++iterBuf )
              {
                CMP::RecvBuffer& recv_buffer = iterBuf->get_message_buffer();
                recv_buffer.irecv();
                //bool flag = recv_buffer.test();
                //printf("reception ID OLD     %d  \n", nstep ); fflush(0);
              }
            for (E_Int ircv = 1; ircv < nbcomID_U +1; ++ircv)
              {
                 E_Int source = param_int_tc[ pt_debID_U + ircv];
                 pt_rcv_queue->emplace_back( source , etiquette);
                 CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
                 recv_buffer.irecv();
                 //bool flag = recv_buffer.test();
                //printf("reception ID Unsteady source  %d %d \n", source, nstep ); fflush(0);
              }
        }

      #ifdef TimeShow
      #ifdef _OPENMP
       E_Float time_out = omp_get_wtime();
       ipt_timecount[0] = ipt_timecount[0] + time_out -time_in;
      #endif
      #endif

      //MPI_Barrier(MPI_COMM_WORLD);

      E_Int nb_send_buffer = 0;
      for (E_Int ip2p = 1; ip2p < Nbp2p_send +1; ++ip2p)
      {
        E_Int ech  = param_int_tc[ip2p + shift_graph];
        dest       = param_int_tc[ech];

        //printf("avt setInterpTransfersInter: dest rank %d %d \n", dest, rank);fflush(0);

        if (dest != rank && param_real_tc != NULL)  // Inter Process ibc
        {
          nb_send_buffer += 1;
          TypeTransfert = 1;
          #ifdef _MPI
          send_mpi=1;
          
          if     (Nopass==0 ){ 
                              K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
                                          param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                          pair_of_queue_pass1, etiquette,  ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer);
                             }
 
          else if(Nopass==1 ){ 
                              K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
                                          param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                          pair_of_queue_pass2, etiquette,  ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer);
                             }
          else if(Nopass==2 ){ 
                              K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
                                          param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                          pair_of_queue_pass3, etiquette,  ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer);
                             }
          else if(Nopass==3 ){ 
                              K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
                                          param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                          pair_of_queue_pass4, etiquette,  ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer);
                             }

          // printf("Transfert inter: dest  %d , nstep,  %d , Nopass %d \n", dest, nstep, Nopass ); fflush(0);
          //K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
          //                                 param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
          //                                 pair_of_queue_loc, etiquette,  ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer);
          #endif
        }
      }//loop comm p2p
    }  //mpi
    // _MPI
    #endif

    //comm local pour recouvrememnt transfert IBC
    for (E_Int ip2p = 1; ip2p < Nbp2p_send +1; ++ip2p)
    {
      E_Int ech  = param_int_tc[ip2p+shift_graph];
      dest       = param_int_tc[ech];
      if (dest == rank && param_real_tc != NULL )  // Intra Process
      { TypeTransfert = 1;
        //printf("Transfert intra: dest  %d , nstep,  %d , Nopass %d  \n", dest, nstep, Nopass ); fflush(0);
        K_FASTC::setInterpTransfersIntra(iptro_tmp    , vartype    , param_int_tc, param_real_tc,
                                         param_int    , param_real , ipt_omp     , linelets_int, linelets_real,  it_target, nidom, ip2p, 
                                         ipt_timecount, nstep      , nitmax      , rk           , exploc       , numpassage); 
      }
    } //loop ip2p
/*
*/

    #ifdef TimeShow
    #ifdef _OPENMP
      time_in = omp_get_wtime();
    #endif
    #endif

    #ifdef _MPI
     if (mpi )
    {
      //comm multi processus: wait + remplissage point cible
      
      E_Int nbcomID = nbcomID_S + nbcomID_U;  // nbr dechange PaP  en reception
      //printf("get ID: nbP2P: %d , nstep: %d ,npass: %d \n", nbcomID,  nstep, Nopass ); fflush(0);
      if (Nopass==0)
          {  if(nbcomID  !=0){ K_FASTC::getTransfersInter(nbcomID, iptro_tmp, param_int, param_real, param_int_tc ,pair_of_queue_pass1, ipt_timecount);}
             if(send_mpi ==1){   SendQueue* pt_snd_queue =  pair_of_queue_pass1->second; pt_snd_queue->waitAll(); pt_snd_queue->clear(); }
          }
      else if (Nopass==1)
          {  if(nbcomID  !=0){ K_FASTC::getTransfersInter(nbcomID, iptro_tmp, param_int, param_real, param_int_tc ,pair_of_queue_pass2, ipt_timecount);}
             if(send_mpi ==1){   SendQueue* pt_snd_queue =  pair_of_queue_pass2->second; pt_snd_queue->waitAll(); pt_snd_queue->clear(); }
          }
      else if (Nopass==2)
          {  if(nbcomID  !=0){ K_FASTC::getTransfersInter(nbcomID, iptro_tmp, param_int, param_real, param_int_tc ,pair_of_queue_pass3, ipt_timecount);}
             if(send_mpi ==1){   SendQueue* pt_snd_queue =  pair_of_queue_pass3->second; pt_snd_queue->waitAll(); pt_snd_queue->clear(); }
          }
      else if (Nopass==3)
          {  if(nbcomID  !=0){ K_FASTC::getTransfersInter(nbcomID, iptro_tmp, param_int, param_real, param_int_tc ,pair_of_queue_pass4, ipt_timecount);}
             if(send_mpi ==1){   SendQueue* pt_snd_queue =  pair_of_queue_pass4->second; pt_snd_queue->waitAll(); pt_snd_queue->clear(); }
          }
      //printf("apres get ID: nbP2P: %d , nstep: %d ,npass: %d \n", nbcomID,  nstep, Nopass ); fflush(0);

      //K_FASTC::getTransfersInter(nbcomID, iptro_tmp, param_int, param_real, param_int_tc , pair_of_queue_loc);

      #ifdef TimeShow
      #ifdef _OPENMP
       E_Float time_out = omp_get_wtime();
       ipt_timecount[4] = ipt_timecount[4] + time_out -time_in;
       time_in= omp_get_wtime();
      #endif
      #endif
    }

/*
    //comm local pour recouvrememnt transfert IBC
    for (E_Int ip2p = 1; ip2p < Nbp2p_send +1; ++ip2p)
    {
      E_Int ech  = param_int_tc[ip2p+shift_graph];
      dest       = param_int_tc[ech];
      if (dest == rank)  // Intra Process
      { TypeTransfert = 1;
        //printf("Transfert intra: dest  %d , nstep,  %d , Nopass %d time alloc %f \n", dest, nstep, Nopass, ipt_timecount[5] ); fflush(0);
        K_FASTC::setInterpTransfersIntra(iptro_tmp    , vartype    , param_int_tc, param_real_tc,
                                         param_int    , param_real , ipt_omp     , linelets_int, linelets_real,  it_target, nidom, ip2p, 
                                         ipt_timecount, nstep      , nitmax      , rk           , exploc       , numpassage); 
      }
    } //loop ip2p

    #ifdef TimeShow
    #ifdef _OPENMP
      time_in = omp_get_wtime();
    #endif
    #endif
*/

    
       #ifdef TimeShow
       #ifdef _OPENMP
        outputfile << "Time in getTransfersInter "     << ipt_timecount[4] <<  " nstepn " << nstep << " npass" << Nopass << std::endl;
        outputfile << "Time InterpTransfert (Intra)  " << ipt_timecount[1] << std::endl;
        outputfile << "Time AllocTransfert (Intra)   " << ipt_timecount[5] << std::endl;
        outputfile << "Time in MPI send_buffer, irecv "<< ipt_timecount[0] << std::endl;
        outputfile << "Time InterpTransfert (Inter)  " << ipt_timecount[2] << std::endl;
        outputfile << "Time Wait Transfert (Inter)   " << ipt_timecount[6] << std::endl;
        outputfile << "Nb com. p2p " << Nbp2p_send << std::endl;
        outputfile << std::endl << std::endl;
        outputfile.close();

        time_in = omp_get_wtime();
       #endif
       #endif
    #endif

  } //if  param_int_tc != Null


}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
void K_FASTC::setInterpTransfersIntra(

    E_Float**& ipt_ro, E_Int& varType, E_Int*& param_int_tc,
    E_Float*& param_real_tc, E_Int**& param_int, E_Float**& param_real,  E_Int*& ipt_omp, E_Int*& linelets_int, E_Float*& linelets_real,
    E_Int& it_target, E_Int& nidom, E_Int& NoTransfert,
    E_Float*& ipt_timecount, E_Int& nstep, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& num_passage)
{

#ifdef TimeShow
#ifdef _OPENMP
  E_Float time_in = omp_get_wtime();
#endif
#endif

  // E_Int NoTransfert   = 1; // ONLY INTRA
  E_Int nvars;
  if      (varType <= 3 && varType >= 1) nvars = 5;
  else if (varType == 4)  nvars = param_int[0][NEQ_LBM] +5;     // LBM transfer, 19 or 27 Qs and 5 macros (32 max in total)
  else if (varType == 41) nvars = param_int[0][NEQ_LBM] +5 + 12;
  else if (varType == 5 ) nvars = param_int[0][NEQ_LBM] +5 +  6; // LBM Overset ou hybride : 19 or 27 Q, 5 macro and 6 gradients
  else                    nvars = 6;

  E_Int* ipt_cnd = NULL;  // ONLY FOR STRUCTURED
  
  E_Int sizecomID  = param_int_tc[2];
  E_Int shift_graph = sizecomID + 2;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech           = param_int_tc[NoTransfert + shift_graph];
  E_Int nrac          = param_int_tc[ech + 1];  // nb total de raccord
  E_Int nrac_inst     = param_int_tc[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel     = param_int_tc[ech + 3];  // nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady   = nrac - nrac_inst;        // nb total de raccord stationnaire

  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;

  if (nrac_inst > 0) pass_inst_fin=2;

  //on optimise les transfert pour implicit local
  E_Int impli_local[nidom];
  if (ipt_omp == NULL) // Merde fastP a gerer
    { for (E_Int nd = 0; nd < nidom; nd++) {impli_local[nd]=1;}  }
  else{ 
        E_Int nbtask = ipt_omp[nstep-1]; 
        E_Int ptiter = ipt_omp[nssiter+ nstep-1];

        for (E_Int nd = 0; nd < nidom; nd++) {impli_local[nd]=0;}
        for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
          E_Int pttask = ptiter + ntask*(6+threadmax_sdm*7);
          E_Int nd = ipt_omp[ pttask ];
          impli_local[nd]=1;
        }
      }

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(  size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation]; // Pour l explicite local

  E_Int ntab_int     =18;
  E_Float cutoff_coef=1.e-12;

  //E_Int rank=0;
  //MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  // on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx = 0; E_Int ibcTypeMax=0;
  for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
  {
    E_Int irac_deb = 0;
    E_Int irac_fin = nrac_steady;
    if(pass_inst == 1){ irac_deb = param_int_tc[ ech + 4 + it_target ]; irac_fin = param_int_tc[ ech + 4 + it_target + timelevel ]; }

   for (E_Int irac = irac_deb; irac < irac_fin; irac++)
    {
      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;

      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      //  si incompatibilite pass et typeTransfert, on skippe le raccord
      E_Int ibcType = param_int_tc[shift_rac + nrac * 3];

      if (ibcType > ibcTypeMax){ ibcTypeMax= ibcType;}
      E_Int ibc = 1;
      if (ibcType < 0) ibc = 0;


      // Si on est en explicit local, on va autoriser les transferts entre certaines zones en fonction de la ss-ite courante
      if(exploc == 1)
	{
	  E_Int debut_rac = ech + 4 + timelevel*2 + nrac*ntab_int + 27*irac;
	  E_Int levelD = param_int_tc[debut_rac + 25];
	  E_Int levelR = param_int_tc[debut_rac + 24];
	  E_Int cyclD  = nssiter/levelD;

	  // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	  if (levelD > levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==cyclD-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1)) { autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	  // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	  else if (levelD < levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)
		{ autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	  // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	  else if (levelD == levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1)
		{ autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	  // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
	  else if (levelD == levelR and num_passage == 2)
	    {
	      if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1) { autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	}
      // Sinon, on autorise les transferts  si la zone donneuse a ete modifiee a l'iteration nstep
      else {
             E_Int NoD      =  param_int_tc[ shift_rac + nrac*5     ];
             if (impli_local[NoD]==1) autorisation_transferts[pass_inst][irac_auto]=1;
           }

      if (autorisation_transferts[pass_inst][irac_auto]==1)
	   { 
	      E_Int nbRcvPts = param_int_tc[shift_rac + nrac*10 ];
    
	      if (nbRcvPts > nbRcvPts_mx) nbRcvPts_mx = nbRcvPts;

           }// autorisation transfert

    }
  }

  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 17 * threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

  E_Float** RcvFields = new E_Float*[ nvars*threadmax_sdm];
  E_Float** DnrFields = new E_Float*[ nvars*threadmax_sdm];

#ifdef TimeShow
#ifdef _OPENMP
  E_Float time_out = omp_get_wtime();
  ipt_timecount[5] = ipt_timecount[5] + time_out - time_in;
  time_in = omp_get_wtime();
#endif
#endif

//# pragma omp parallel default(shared)  num_threads(1)
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int ithread           = omp_get_thread_num() + 1;
    E_Int Nbre_thread_actif = omp_get_num_threads();  // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc /*, nocf*/, indCoef, noi, sizecoefs,
        /*Nbchunk,*/ imd, jmd, imdjmd;

    E_Float** vectOfRcvFields = RcvFields + nvars*(ithread-1);
    E_Float** vectOfDnrFields = DnrFields +  nvars*(ithread-1);

    // 1ere pass_typ: IBC
    // 2eme pass_typ: transfert
    //

    E_Int count_racIBC = 0;

    // 1ere pass_inst: les raccord fixe
    // 2eme pass_inst: les raccord instationnaire
    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
      {
        //  printf("pass %d %d %d %d \n", rank ,ipass_typ, pass_inst, it_target ); 
        // printf("ipass_inst = %d, level= %d \n",  ipass_inst, nrac_inst_level );

        E_Int irac_deb = 0;
        E_Int irac_fin = nrac_steady;
        if (pass_inst == 1) {
          irac_deb = param_int_tc[ech + 4 + it_target];
          irac_fin = param_int_tc[ech + 4 + it_target + timelevel];
        }

        // printf("iracdeb=  %d, iracfin= %d \n", irac_deb, irac_fin  );
        for (E_Int irac = irac_deb; irac < irac_fin; irac++) {

	 E_Int irac_auto= irac-irac_deb;
	 if (autorisation_transferts[pass_inst][irac_auto]==1)
	  {

          E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
          //printf("ipass_typ = %d, pass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,pass_inst,irac , ithread );
          // ipass_typ,ipass_inst,irac , ithread );
          E_Int ibcType = param_int_tc[shift_rac + nrac * 3];
          E_Int ibc = 1;
          if (ibcType < 0) ibc = 0;

          E_Int NoD       = param_int_tc[shift_rac + nrac * 5  ];
          E_Int NoR       = param_int_tc[shift_rac + nrac * 11 ];
          E_Int nvars_loc = param_int_tc[shift_rac + nrac * 13 ];  // neq fonction raccord rans/LES
          E_Int rotation  = param_int_tc[shift_rac + nrac * 14 ];  // flag pour periodicite azymuthal

          // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
          E_Int solver_D=2; E_Int solver_R=2;
          if (nvars_loc == 11) {solver_R =4;}
          if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
          if (nvars_loc == 19) {solver_D =4; solver_R=4;}

          E_Int overset  =  param_int[NoD][LBM_OVERSET];        //flag pour overset en LBM
          if      (nvars_loc== param_int[NoD][NEQ_LBM] && overset==0) nvars_loc = nvars_loc + 5;
          else if (nvars_loc== param_int[NoD][NEQ_LBM] && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;
 
          //printf("nvar loc %d , solver_RD= %d %d \n", nvars_loc, solver_R, solver_R);

          E_Int meshtype = param_int[ NoD ][ MESHTYPE ] ;

          E_Int cnNfldD = 0;
          E_Int* ptrcnd = NULL;

          if (nvars_loc == 5 || nvars_loc == 6) // Transferts NS classiques ou LBM -> NS
          {
           for (E_Int eq = 0; eq < nvars_loc; eq++)
           {
             vectOfRcvFields[eq] = ipt_ro[NoR] + eq * param_int[NoR][ NDIMDX];
             vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX];
           }
          }
          else if (nvars_loc == param_int[NoD][NEQ_LBM] + 5) // Transferts LBM classiques
          {
            // On commence par copier les 5 variables macros
            for (E_Int eq = 0; eq < 5; eq++)
            {
              vectOfRcvFields[eq] = ipt_ro[NoR] + eq * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX];
            }
            // Puis on copie les fonctions de distribution
            for (E_Int eq = 5; eq < nvars_loc; eq++)
            {
              vectOfRcvFields[eq] = ipt_ro[NoR + nidom] + (eq-5) * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD + nidom] + (eq-5) * param_int[NoD][ NDIMDX];
            }
          }
          else if (nvars_loc == 11 ) // //Transfert NS -> LBM    
          {
            // On commence par copier les 5 variables macros
            for (E_Int eq = 0; eq < 5; eq++)
            {
              vectOfRcvFields[eq] = ipt_ro[NoR] + eq * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX];
            }
            // Puis on copie les gradients
            for (E_Int eq = 5; eq < nvars_loc; eq++)
            {
              vectOfRcvFields[eq] = ipt_ro[NoR + nidom] + (eq-5) * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD + nidom] + (eq-5) * param_int[NoD][ NDIMDX];
            }
          }
          else if (nvars_loc == param_int[NoD][NEQ_LBM] + 17 ) // //Transfert LBM  overset   
          {
            // On commence par copier les 5 variables macros
            for (E_Int eq = 0; eq < 5; eq++)
            {
              vectOfRcvFields[eq] = ipt_ro[NoR] + eq * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX];
            }
            // Puis on copie les gradients
            for (E_Int eq = 5; eq < 11; eq++) //S
            {
              vectOfRcvFields[eq] = ipt_ro[NoR + nidom*2] + (eq-5) * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD + nidom*2] + (eq-5) * param_int[NoD][ NDIMDX];
            }
            for (E_Int eq =11; eq < 17; eq++) //psiG
            {
              vectOfRcvFields[eq] = ipt_ro[NoR + nidom*3] + (eq-11) * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD + nidom*3] + (eq-11) * param_int[NoD][ NDIMDX];
            }
            for (E_Int eq =17; eq < nvars_loc; eq++) //Q
            {
              vectOfRcvFields[eq] = ipt_ro[NoR + nidom] + (eq-17) * param_int[NoR][ NDIMDX];
              vectOfDnrFields[eq] = ipt_ro[NoD + nidom] + (eq-17) * param_int[NoD][ NDIMDX];
            }
          }

          imd = param_int[ NoD ][NIJK  ];
          jmd = param_int[ NoD ][NIJK+1];
          imdjmd = imd * jmd;

          ////
          //  Interpolation parallele
          ////
          ////

          E_Int nbRcvPts = param_int_tc[shift_rac + nrac*10];
          // E_Int nbDonPts = param_int_tc[ shift_rac                ];

          E_Int pos;
          pos = param_int_tc[shift_rac + nrac*7 ]; E_Int* ntype    = param_int_tc  + pos;
          pos = pos + 1 + ntype[0];                E_Int* types    = param_int_tc  + pos;
          pos = param_int_tc[shift_rac + nrac*6 ]; E_Int* donorPts = param_int_tc  + pos;
          pos = param_int_tc[shift_rac + nrac*12]; E_Int* rcvPts   = param_int_tc  + pos;  // donor et receveur inverser car storage donor
          pos = param_int_tc[shift_rac + nrac*8 ]; E_Float* ptrCoefs = param_real_tc + pos;

          E_Int nbInterpD = param_int_tc[shift_rac + nrac];
          E_Int nbFluCons = param_int_tc[shift_rac + nrac*16];

          // sommation flux fin et stokage dans flux grossier pour nearmatch conservatif
          for (E_Int nbflu = 0; nbflu < nbFluCons; nbflu++)
            {
             if(nvars_loc==5)
              {
               #include "FastC/Com/flux_conservatif_5eq.cpp"
              }
             else if(nvars_loc==6)
              {
               #include "FastC/Com/flux_conservatif_6eq.cpp"
              }
            }
          E_Int ideb = 0;
          E_Int ifin = 0;
          E_Int shiftCoef = 0;
          E_Int shiftDonor = 0;

          for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++) {
            type = types[ifin];

            SIZECF(type, meshtype, sizecoefs);
            ifin = ifin + ntype[1 + ndtyp];


            // // *      New school: meilleur equilibrage, mais gestion looop
            // dynamique rame...
            // // *
            //         E_Int size_bc =  ifin-ideb;
            //         E_Int size_min=   16      ;
            //         //E_Int chunk = size_bc/Nbre_thread_actif;
            //         //if     (chunk < size_min && size_bc >= size_min) { chunk =
            // size_min;}
            //         //else if(chunk < size_min && size_bc <  size_min) { chunk =
            // size_bc ;} E_Int chunk = size_min; if(size_bc <  size_min) { chunk =
            // size_bc ;}

            //         if      ( type == 0 ||  chunk <= 0) { Nbchunk = 1; } else if
            // ( chunk > 0)                { Nbchunk = size_bc/chunk;}

            //         chunk = size_bc/Nbchunk;

            //         E_Int r = size_bc - chunk*Nbchunk;

            //         #pragma omp for nowait schedule(dynamic,1)
            //         for (E_Int nd = 0; nd < Nbchunk; nd++)
            //         {

            E_Int pt_deb, pt_fin;

            /// oldschool
            // Calcul du nombre de champs a traiter par chaque thread
            E_Int size_bc = ifin - ideb;
            E_Int chunk = size_bc / Nbre_thread_actif;
            E_Int r = size_bc - chunk * Nbre_thread_actif;
            // pts traitees par thread
            if (ithread <= r) {
              pt_deb = ideb + (ithread - 1) * (chunk + 1);
              pt_fin = pt_deb + (chunk + 1);
            } else {
              pt_deb = ideb + (chunk + 1) * r + (ithread - r - 1) * chunk;
              pt_fin = pt_deb + chunk;
            }

            // Si type 0, calcul sequentiel
            if (type == 0) {
              if (ithread == 1) {
                pt_deb = ideb;
                pt_fin = ifin;
              } else {
                pt_deb = ideb;
                pt_fin = ideb;
              }
            }

            /// newschool suite
            //        if (nd  <  r) { pt_deb = ideb + nd*(chunk+1)
            //        ; pt_fin = pt_deb + (chunk+1); } else          { pt_deb =
            //        ideb +    (chunk+1)*r+(nd-r)*chunk; pt_fin = pt_deb +
            //        chunk;    }

            //if(NoR == 36 || NoR == 36) printf(" irac= %d, NoR= %d, nvar=  %d, NoD= %d, Rans=  %d, rot=
            // %d, fin= %d, ithread= %d \n", irac, NoR, nvars_loc, NoD,
            // ipass_inst ,rotation, pt_fin , ithread ); 
            // if(ithread <=8 &&
            // NoD==83 )  printf(" shift %d  %d %d %d  %d %d %d  %d \n", irac,
            // NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );
            // if(ithread <=8 && NoR==114 )  printf(" new   %d  %d %d %d  %d %d
            // %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  ,
            // type, ithread );

                      noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
                      indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
                      E_Int shiftv =0;
                      if     (nvars_loc==5 || (ibc==1 && solver_R==4) )
                      {
            #           include "TRANSFERT/commonInterpTransfers_reorder_5eq.h"
                      }
                      else if(nvars_loc==6)
                      {
            #           include "TRANSFERT/commonInterpTransfers_reorder_6eq.h"
                      }
                      else
                      {
	    #           include "TRANSFERT/commonInterpTransfers_reorder_neq.h"
                      }

                      // COUPLAGE NS-LBM: changement d'unite
                      if (solver_D==4 && solver_R<4)
                      {
                         // Transfert LBM vers NS: repasse dans unites SI
#                        include "TRANSFERT/includeTransfers_dimLBMtoNS.h"
                      }
                      else if (solver_D<4 && solver_R==4)
                      {
                         // Transfert NS vers LBM : adimensionnement
#                        include "TRANSFERT/includeTransfers_dimNStoLBM.h"
                      }


                      // Prise en compte de la periodicite par rotation
                      if (rotation == 1)
                      {
                       E_Float* angle = ptrCoefs + nbInterpD;
            #          include "TRANSFERT/includeTransfers_rotation.h"
                      }
                      // ibc
                      if (ibc == 1)
			{
                          E_Float* linelets    = NULL;
                          E_Int* indexlinelets = NULL;
                          E_Int nbptslinelets  = 0;
                          E_Float* xPC     = ptrCoefs + nbInterpD;
                          E_Float* xPI     = ptrCoefs + nbInterpD + 3 * nbRcvPts;
                          E_Float* xPW     = ptrCoefs + nbInterpD + 6 * nbRcvPts;
                          E_Float* densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;
                          if (linelets_int != NULL )
                              {
                                nbptslinelets        = linelets_int[0];
                                E_Int addrlinelets   = linelets_int[count_racIBC + 3 ];
                                linelets             = linelets_real + addrlinelets;
                                indexlinelets        = linelets_int + linelets_int[1]+1 + 3;
                                count_racIBC         = count_racIBC + 1;
                              }
			  K_FASTC::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                                 xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                                                 xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                                                 xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                                                 densPtr, 
                                                                 ipt_tmp, size, nvars,
                                                                 param_real[ NoR ],
                                                                 vectOfDnrFields, vectOfRcvFields,
                                                                 nbptslinelets, linelets, indexlinelets);
                      }//ibc

                      //*
                      //        } //chunk
                      //*/
                      ideb       = ideb + ntype[1 + ndtyp];
                      shiftCoef  = shiftCoef   +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
                      shiftDonor= shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif

                   }// type
	          } //autorisation transfert
                }//irac
               }//pass_inst
  }    // omp

  delete[] ipt_cnd; delete [] RcvFields;  delete [] DnrFields;

#ifdef TimeShow
#ifdef _OPENMP
  time_out = omp_get_wtime();
  ipt_timecount[1] = ipt_timecount[1] + time_out - time_in;
  time_in = omp_get_wtime();
#endif
#endif

  // return varType;
}

#ifdef _MPI
//=============================================================================
// Transfert de champs sous forme de numpy CMP lib
// From zone
// Retourne une liste de numpy directement des champs interpoles
// in place + from zone + tc compact
//=============================================================================
void K_FASTC::setInterpTransfersInter(
    E_Float**& ipt_ro   , E_Int& varType   , E_Int*& param_int_tc, E_Float*& param_real_tc,
    E_Int**& param_int  , E_Float**& param_real, E_Int*& ipt_omp, E_Int*& linelets_int    , E_Float*& linelets_real, 
    E_Int& TypeTransfert, E_Int& it_target, E_Int& nidom, E_Int& NoTransfert,
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue_loc, E_Int& etiquette,
    E_Float*& ipt_timecount                          ,
    E_Int& nstep, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& num_passage, E_Int& nb_send_buffer)

{
#ifdef TimeShow
#ifdef _OPENMP
    E_Float time_in = omp_get_wtime();
#endif
#endif

  E_Int nvars;
  if      (varType <= 3 && varType >= 1) nvars = 5;
  else if (varType == 4)  nvars = param_int[0][NEQ_LBM] +5;     // LBM transfer, 19 or 27 Qs and 5 macros (32 max in total)
  else if (varType == 41) nvars = param_int[0][NEQ_LBM] +5 +12; // LBM Overset ou hybride : 19 or 27 Q, 5 macro and 6 gradients
  else if (varType == 5 ) nvars = param_int[0][NEQ_LBM] +5 + 6; // LBM Overset ou hybride : 19 or 27 Q, 5 macro and 6 gradients
  else                    nvars = 6;

  E_Int sizecomID = param_int_tc[2];
  E_Int shift_graph = sizecomID + 2;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech       = param_int_tc[NoTransfert + shift_graph];
  E_Int dest      = param_int_tc[ech    ];  // processus destination
  E_Int nrac      = param_int_tc[ech + 1];  // nb total de raccord
  E_Int nrac_inst = param_int_tc[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel = param_int_tc[ech + 3];  // nb de pas de temps stocker pour
                                             // chaque raccord instationnaire

  E_Int nrac_steady = nrac - nrac_inst;  // nb total de raccord stationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;
  E_Int nrac_inst_level = 0;
  if (nrac_inst > 0) {
   pass_inst_fin=2;
   nrac_inst_level = param_int_tc[ech + 4 + it_target + timelevel] - param_int_tc[ech + 4 + it_target] + 1; 
  } 
  
  //printf("send %d , nrac= %d , nrac_inst= %d , timelevelNb= %d \n", dest, nrac, nrac_inst, nrac_inst_level);

  // on dimension tableau travail pour IBC et pour transfert
  // E_Int nrac_inst_level = param_int_tc[ech + 4 + it_target + timelevel] -
  //                         param_int_tc[ech + 4 + it_target] + 1;

  std::vector<E_Float*> frp(nrac_steady + nrac_inst_level);

  int rank;
#ifdef _MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#endif
  SendQueue* pt_snd_queue = pair_of_queue_loc->second;
  if (pt_snd_queue->size() < (size_t)nb_send_buffer)
  {
      pt_snd_queue->emplace_back(dest, etiquette);

      //printf("size sendQeueNEW  %d %d dest= %d , etiquetee= %d \n",pt_snd_queue->size(), nb_send_buffer, dest, etiquette); 
  }

  // Preparation du buffer d'envoi :
  CMP::SendBuffer& send_buffer = *(*pt_snd_queue)[nb_send_buffer-1].message_buffer;//back_message_buffer();
  
  // A partir d'ici pour allouer les tableaux a remplir
  std::vector<CMP::SendBuffer::PackedData*> pck_data;

  //on optimise les transfert pour implicit local
  E_Int impli_local[nidom];
  if (ipt_omp == NULL) // Merde fastP a gerer
    { for (E_Int nd = 0; nd < nidom; nd++) {impli_local[nd]=1;}  }
  else{ 
        E_Int nbtask = ipt_omp[nstep-1]; 
        E_Int ptiter = ipt_omp[nssiter+ nstep-1];

        for (E_Int nd = 0; nd < nidom; nd++) {impli_local[nd]=0;}
        for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
          E_Int pttask = ptiter + ntask*(6+threadmax_sdm*7);
          E_Int nd = ipt_omp[ pttask ];
          impli_local[nd]=1;
        }
      }

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(  size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  //E_Int autorisation_transferts[500][500];
  // Test IF Has Data to send and count rac:
  // Nb rac.
  bool has_data_to_send = false;
  E_Int count_rac = 0;
  E_Int nbRcvPts_mx = 0;
  E_Int ibcTypeMax  = 0;
  E_Int ntab_int    =18;
  E_Float cutoff_coef=1.e-12;

  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
        E_Int irac_deb = 0;
        E_Int irac_fin = nrac_steady;
        if ( pass_inst == 1 )
        {
            irac_deb = param_int_tc[ech + 4 + it_target];
            irac_fin = param_int_tc[ech + 4 + it_target + timelevel];
        }

        for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
        {

	  E_Int irac_auto= irac-irac_deb;
	  autorisation_transferts[pass_inst][irac_auto]=0;

	  E_Int shift_rac = ech + 4 + timelevel * 2 + irac;

          //  si incompatibilite pass et typeTransfert, on skippe le raccord
          E_Int ibcType = param_int_tc[shift_rac + nrac*3];

          if (ibcType > ibcTypeMax){ ibcTypeMax= ibcType;}
          E_Int ibc = 1;
	  if (ibcType < 0) ibc = 0;

	  if(exploc == 1)
	    {
	      E_Int debut_rac = ech + 4 + timelevel*2 + nrac*ntab_int + 27*irac;
	      E_Int levelD = param_int_tc[debut_rac + 25];
	      E_Int levelR = param_int_tc[debut_rac + 24];
	      E_Int cyclD  = nssiter/levelD;

	      // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	      if (levelD > levelR and num_passage == 1)
		{
		  if (nstep%cyclD==cyclD-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1)){autorisation_transferts[pass_inst][irac_auto]=1;}
		}
	      // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	      else if (levelD < levelR and num_passage == 1)
		{
		  if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)
                     {autorisation_transferts[pass_inst][irac_auto]=1;}
		}
	      // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	      else if (levelD == levelR and num_passage == 1)
		{
		  if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1) { autorisation_transferts[pass_inst][irac_auto]=1; }

		}
	      // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
	      else if (levelD ==  levelR and num_passage == 2)
		{
		  if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1){autorisation_transferts[pass_inst][irac_auto]=1;}
		}
	    }
           // Sinon, on autorise les transferts  si la zone donneuse a ete modifiee a l'iteration nstep
	   else 
        {
            //E_Int NoD      =  param_int_tc[ shift_rac + nrac*5     ];
            //if (impli_local[NoD]==1) autorisation_transferts[pass_inst][irac_auto]=1;
            autorisation_transferts[pass_inst][irac_auto]=1;
        }
	  
	  if (autorisation_transferts[pass_inst][irac_auto]==1)
	   {
	      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	      E_Int nbRcvPts = param_int_tc[shift_rac + nrac*10];

	      if (nbRcvPts > nbRcvPts_mx) nbRcvPts_mx = nbRcvPts;
	      has_data_to_send = true;
	      count_rac += 1;

	     }// autorisation transfert
        } // irac
    }  // pas inst

if (has_data_to_send)
{
#ifdef TimeShow
#ifdef _OPENMP
 E_Float time_out = omp_get_wtime();
 ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
 time_in  = omp_get_wtime();
#endif
#endif

    send_buffer << count_rac;

    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) {
          irac_deb = param_int_tc[ech + 4 + it_target];
          irac_fin = param_int_tc[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
      {

	E_Int irac_auto= irac-irac_deb;
	if (autorisation_transferts[pass_inst][irac_auto]==1)
	{ 
	    E_Int shift_rac = ech + 4 + timelevel * 2 + irac;

	    E_Int nbRcvPts  = param_int_tc[shift_rac + nrac*10];
	    E_Int nvars_loc = param_int_tc[shift_rac + nrac*13];  // flag raccord rans/LES
	    E_Int Nozone    = param_int_tc[shift_rac + nrac*11];
            E_Int nbFluCons = param_int_tc[shift_rac + nrac*16];

            // on determine un No zone pipeau pour skipper remplissage inutile en implicit local
            E_Int NoD       = param_int_tc[shift_rac + nrac*5];
            E_Int Nozone_loc   = Nozone; 
            E_Int nbRcvPts_loc = nbRcvPts;
            if (impli_local[NoD] == 0) {Nozone_loc = -999; nbRcvPts_loc=1;}

            // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
            E_Int solver_D=2; E_Int solver_R=2;
            if (nvars_loc == 11) {solver_R =4;}
            if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
            if (nvars_loc == 19) {solver_D =4; solver_R=4;}

            E_Int overset  =  param_int[NoD][LBM_OVERSET];        //flag pour overset en LBM
            if      (nvars_loc== param_int[NoD][NEQ_LBM] && overset==0) nvars_loc = nvars_loc + 5;
            else if (nvars_loc== param_int[NoD][NEQ_LBM] && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;

            //if(Nozone==5 and count_rac==0 ) { printf("Nozone_loc %d \n", Nozone_loc);}

	    send_buffer << Nozone_loc;
           
            E_Int  nbFlu_nvar = nbFluCons +100*nvars_loc;
	    send_buffer << nbFlu_nvar;

            //dimensionnement tableau real pour envoi (interp + flux)
            E_Int sz_real=nbRcvPts_loc;
            for (E_Int nflu=0; nflu < nbFluCons; nflu++)
              { E_Int   pos = param_int_tc[shift_rac + nrac*17];
                E_Int* fluD = param_int_tc + pos +7*nflu;
                E_Int idir    = fluD[0];
                E_Int sizefluD= fluD[1];
                E_Int ratio_i = fluD[4];
                E_Int ratio_j = fluD[5];
                E_Int ratio_k = fluD[6];
                E_Int scale; //rapport nombre de face fin/grossier
                if      (idir <=2) { scale = ratio_j*ratio_k;}
                else if (idir <=4) { scale = ratio_i*ratio_k;}
                else               { scale = ratio_i*ratio_j;}
                sz_real += sizefluD/scale;
              }
            
	    //pck_data.push_back(&send_buffer.push_inplace_array(nvars_loc * nbRcvPts_loc * sizeof(E_Float) ));    
	    pck_data.push_back(&send_buffer.push_inplace_array(nvars_loc * sz_real * sizeof(E_Float) ));    
	    E_Int PtlistDonor  = param_int_tc[shift_rac + nrac*12];
	    E_Int* ipt_listRcv = param_int_tc + PtlistDonor;

	    send_buffer << CMP::vector_view<E_Int>(ipt_listRcv, nbRcvPts_loc);
 
            if(nbFluCons !=0)
             { E_Int   pos = param_int_tc[shift_rac + nrac*17];
               E_Int* fluD = param_int_tc + pos;
	       send_buffer << CMP::vector_view<E_Int>(fluD, 7*nbFluCons);
             }
            //printf("size rac:  nbRcvPts_loc*nvars_loc: %d , count_rac: %d , irac: %d , passinst: %d \n", sz_real*nvars_loc, count_rac, irac, pass_inst);

	 } // autorisation transfert
      } // irac
    } // pass_inst


    send_buffer.finalize_and_copy();

    count_rac = 0;
    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) {
          irac_deb = param_int_tc[ech + 4 + it_target];
          irac_fin = param_int_tc[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
      {
	E_Int irac_auto= irac-irac_deb;
	if (autorisation_transferts[pass_inst][irac_auto]==1)
	 {
	    frp[count_rac] = pck_data[count_rac]->data<E_Float>();
	    count_rac += 1;
	 }// autorisation transfert
      } // irac
    } // pass_inst
  } // if has data to send


#ifdef TimeShow
#ifdef _OPENMP
    E_Float time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
    time_in  = omp_get_wtime();
#endif
#endif
    
  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 17 * threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

  // // tableau temporaire pour utiliser la routine commune
  // K_CONNECTOR::setIBCTransfersCommon
  FldArrayI rcvPtsI(nbRcvPts_mx);
  E_Int* rcvPts = rcvPtsI.begin();

  E_Float** RcvFields = new E_Float*[ nvars*threadmax_sdm];
  E_Float** DnrFields = new E_Float*[ nvars*threadmax_sdm];

//# pragma omp parallel default(shared)  num_threads(1)
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int ithread = omp_get_thread_num() + 1;
    E_Int Nbre_thread_actif =
        omp_get_num_threads();  // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int type;
    E_Int indD0, indD, i, j, k, ncfLoc, indCoef, noi, sizecoefs, imd, jmd, imdjmd;

    E_Float** vectOfRcvFields = RcvFields + nvars*(ithread-1);
    E_Float** vectOfDnrFields = DnrFields + nvars*(ithread-1);

    // 1ere pass: IBC
    // 2eme pass: transfert
    //

    E_Int count_racIBC = 0;

    // 1ere pass_inst: les raccord fixe
    // 2eme pass_inst: les raccord instationnaire
    E_Int count_rac = 0;

    for ( E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
            {
                E_Int irac_deb = 0;
                E_Int irac_fin = nrac_steady;
                if ( pass_inst == 1 )
                {
                    irac_deb = param_int_tc[ech + 4 + it_target];
                    irac_fin = param_int_tc[ech + 4 + it_target + timelevel];
                }

                for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
                {

                  E_Int shift_rac = ech + 4 + timelevel*2 + irac;
                  E_Int NoD       = param_int_tc[shift_rac + nrac*5 ];
                  E_Int nvars_loc = param_int_tc[shift_rac + nrac*13];  // neq fonction raccord rans/LES


		  E_Int irac_auto = irac-irac_deb;

                  // on determine un envoi pipeau de taille 1*neq pour skipper remplissage inutile en implicit local
		  if ( impli_local[NoD]==0 and autorisation_transferts[pass_inst][irac_auto]==1)
                    {
                      //for (E_Int eq = 0; eq < nvars_loc; eq++) { frp[count_rac][eq] =0; }
                      count_rac += 1;
                    }

                  if (autorisation_transferts[pass_inst][irac_auto]==1 and impli_local[NoD]==1)
		   {

                     E_Int ibcType = param_int_tc[shift_rac + nrac*3];

                     E_Int ibc = 1;
	             if (ibcType < 0) ibc = 0;

                     E_Int nbRcvPts  = param_int_tc[shift_rac + nrac*10];
                     E_Int NoR       = param_int_tc[shift_rac + nrac*11];
                     E_Int rotation  = param_int_tc[shift_rac + nrac*14];  // flag pour periodicite azymuthal

                     // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
                     E_Int solver_D=2; E_Int solver_R=2;
                     if (nvars_loc == 11) {solver_R =4;}
                     if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
                     if (nvars_loc == 19) {solver_D =4; solver_R=4;}

                     E_Int overset  =  param_int[NoD][LBM_OVERSET];        //flag pour overset en LBM
                     if      (nvars_loc== param_int[NoD][NEQ_LBM] && overset==0) nvars_loc = nvars_loc + 5;
                     else if (nvars_loc== param_int[NoD][NEQ_LBM] && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;

                     E_Int meshtype = 1;  // ONLY FOR STRUCTURE ipt_ndimdxD[NoD + nidom*6];
                     E_Int cnNfldD  = 0;
                     E_Int* ptrcnd  = NULL;


                     // printf("navr_loc %d %d %d \n", nvars_loc, nvars, Rans);

                     if (nvars_loc == 5 || nvars_loc == 6) // Transferts NS classiques ou LBM -> NS
                     {
                      for (E_Int eq = 0; eq < nvars_loc; eq++) { vectOfDnrFields[eq] = ipt_ro[NoD]    + eq * param_int[NoD][ NDIMDX]; }
                     }
                     else if (nvars_loc == param_int[NoD][NEQ_LBM] + 5) // Transferts LBM classiques
                     {
                       // On commence par copier les 5 variables macros
                       for (E_Int eq = 0; eq < 5; eq++)         { vectOfDnrFields[eq] = ipt_ro[NoD]         +     eq * param_int[NoD][ NDIMDX]; }
                       // Puis on copie les fonctions de distribution
                       for (E_Int eq = 5; eq < nvars_loc; eq++) { vectOfDnrFields[eq] = ipt_ro[NoD + nidom] + (eq-5) * param_int[NoD][ NDIMDX]; }
                     }
                     else if (nvars_loc == 11 ) // //Transfert NS -> LBM    
                     {
                       // On commence par copier les 5 variables macros
                       for (E_Int eq = 0; eq < 5; eq++)         { vectOfDnrFields[eq] = ipt_ro[NoD        ] +     eq * param_int[NoD][ NDIMDX]; }
                       // Puis on copie les gradients
                       for (E_Int eq = 5; eq < nvars_loc; eq++) { vectOfDnrFields[eq] = ipt_ro[NoD + nidom] + (eq-5) * param_int[NoD][ NDIMDX]; }
                     }
                     else if (nvars_loc == param_int[NoD][NEQ_LBM] + 17 ) // //Transfert LBM  overset   
                     {
                       // On commence par copier les 5 variables macros
                       for (E_Int eq = 0; eq < 5; eq++) { vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX]; }
                       // Puis on copie les gradients
                       for (E_Int eq = 5; eq < 11; eq++) //S
                       { vectOfDnrFields[eq] = ipt_ro[NoD + nidom*2] +  (eq-5) * param_int[NoD][ NDIMDX]; }
                       for (E_Int eq =11; eq < 17; eq++) //psiG
                       { vectOfDnrFields[eq] = ipt_ro[NoD + nidom*3] + (eq-11) * param_int[NoD][ NDIMDX]; }
                       for (E_Int eq =17; eq < nvars_loc; eq++) //Q
                       { vectOfDnrFields[eq] = ipt_ro[NoD + nidom  ] + (eq-17) * param_int[NoD][ NDIMDX]; }
                     }

                     for (E_Int eq = 0; eq < nvars_loc; eq++) { vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts; }

                     imd = param_int[ NoD ][ NIJK   ];
                     jmd = param_int[ NoD ][ NIJK+1 ];

                     imdjmd = imd * jmd;
                     ////
                     //  Interpolation parallele
                     ////
                     ////
                     E_Int pos;
                     pos               = param_int_tc[shift_rac + nrac * 7];
                     E_Int* ntype      = param_int_tc + pos;
                     pos               = pos + 1 + ntype[0];
                     E_Int* types      = param_int_tc + pos;
                     pos               = param_int_tc[shift_rac + nrac * 6];
                     E_Int* donorPts   = param_int_tc + pos;
                     pos               = param_int_tc[shift_rac + nrac * 8];
                     E_Float* ptrCoefs = param_real_tc + pos;
  
                     E_Int nbInterpD = param_int_tc[shift_rac + nrac];
                     E_Int nbFluCons = param_int_tc[shift_rac + nrac*16];
                     E_Float* xPC = NULL;
                     E_Float* xPI = NULL;
                     E_Float* xPW = NULL;
                     E_Float* densPtr = NULL;
                     E_Float* linelets   = NULL;
                     E_Int* indexlinelets = NULL;
                     E_Int nbptslinelets = 0;
                     if (ibc == 1) {
                       xPC = ptrCoefs + nbInterpD;
                       xPI = ptrCoefs + nbInterpD + 3 * nbRcvPts;
                       xPW = ptrCoefs + nbInterpD + 6 * nbRcvPts;
                       densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;

                       if (linelets_int != NULL )
                          {
                           nbptslinelets= linelets_int[0];
                           //E_Int addrlinelets   = linelets_int[count_racIBC + 3 ];
                            //E_Float* linelets    = linelets_real + addrlinelets;
                           //E_Int* indexlinelets = linelets_int + linelets_int[1]+1 + 3;
                           count_racIBC = count_racIBC + 1;
                          }
                     }
                     for (E_Int nbflu = 0; nbflu < nbFluCons; nbflu++)
                       { E_Int shift_fluR=0;
                        if(nvars_loc==5)
                        {
                          #include "FastC/Com/flux_conservatif_5eq_D.cpp"
                        }
                        else if(nvars_loc==6)
                         {
                          #include "FastC/Com/flux_conservatif_6eq_D.cpp"
                         }
                       }

                     E_Int ideb = 0;
                     E_Int ifin = 0;
                     E_Int shiftCoef = 0;
                     E_Int shiftDonor = 0;

                     for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++) {
                       type = types[ifin];

                       SIZECF(type, meshtype, sizecoefs);

                       ifin = ifin + ntype[1 + ndtyp];

                       E_Int pt_deb, pt_fin;

                       /// oldschool
                       // Calcul du nombre de champs a traiter par chaque thread
                       E_Int size_bc = ifin - ideb;
                       E_Int chunk = size_bc / Nbre_thread_actif;
                       E_Int r = size_bc - chunk * Nbre_thread_actif;
                       // pts traitees par thread
                       if (ithread <= r) {
                         pt_deb = ideb + (ithread - 1) * (chunk + 1);
                         pt_fin = pt_deb + (chunk + 1);
                       } else {
                         pt_deb = ideb + (chunk + 1) * r + (ithread - r - 1) * chunk;
                         pt_fin = pt_deb + chunk;
                       }

                       // Si type 0, calcul sequentiel
                       if (type == 0) {
                         if (ithread == 1) {
                           pt_deb = ideb;
                           pt_fin = ifin;
                         } else {
                           pt_deb = ideb;
                           pt_fin = ideb;
                         }
                       }
 

                      noi     = shiftDonor;  // compteur sur le tableau d indices donneur
                      indCoef = ( pt_deb - ideb ) * sizecoefs + shiftCoef;
                      E_Int NoR = param_int_tc[shift_rac + nrac*11 ];
                      //if (param_int_tc[ech]==0) printf("No rac= %d , NoR= %d, NoD= %d, Ntype= %d, ptdeb= %d, ptfin= %d, NptD= %d, neq= %d, skip= %d, rank= %d, dest= %d,  thread= %d\n",
                      //irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin , 
                      //param_int_tc[ shift_rac + nrac*10  ], param_int_tc[ shift_rac + nrac*13  ], param_int_tc[ shift_rac + nrac*15  ], 
                      //rank, param_int_tc[ ech  ], ithread );

                      if ( nvars_loc == 5 || (ibc==1 && solver_R==4) ) {
#include "TRANSFERT/commonInterpTransfersD_reorder_5eq.h"
                        } else if ( nvars_loc == 6 ) {
#include "TRANSFERT/commonInterpTransfersD_reorder_6eq.h"
                        } else if ( nvars_loc == 19 ) {
#include "TRANSFERT/commonInterpTransfersD_reorder_19eq.h"
                        } else {
#include "TRANSFERT/commonInterpTransfersD_reorder_neq.h"
                        }

                      // COUPLAGE NS-LBM: changement d'unite
                      if (solver_D==4 && solver_R<4)
                      {
                         // Transfert LBM vers NS: repasse dans unites SI
#                        include "TRANSFERT/includeTransfersD_dimLBMtoNS.h"
                      }
                      else if (solver_D<4 && solver_R==4)
                      {
                         // Transfert NS vers LBM : adimensionnement
#                        include "TRANSFERT/includeTransfersD_dimNStoLBM.h"
                      }

                      // Prise en compte de la periodicite par rotation
                      if ( rotation == 1 ) {
                          E_Float* angle = ptrCoefs + nbInterpD;
#include "TRANSFERT/includeTransfersD_rotation.h"
                          }

                      // ibc
                      if (ibc == 1)
                      {
                        // tableau temporaire pour utiliser la routine commune K_CONNECTOR::setIBCTransfersCommon
                        for ( E_Int noind = pt_deb; noind < pt_fin; noind++ ) rcvPts[noind] = noind;
	  	        K_FASTC::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, 
			                                      xPC, xPC + nbRcvPts, xPC + nbRcvPts * 2, 
							      xPW, xPW + nbRcvPts, xPW + nbRcvPts * 2, 
							      xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
							      densPtr, 
							      ipt_tmp, size, nvars,
							      param_real[ NoD ],
							      vectOfDnrFields, vectOfRcvFields,
							      nbptslinelets, linelets, indexlinelets);
                  
                      }  // ibc
	              //E_Int PtlistDonor  = param_int_tc[shift_rac + nrac*12];
	              //E_Int* ipt_listRcv = param_int_tc + PtlistDonor;

                      //        } //chunk
                      ideb        = ideb + ntype[1 + ndtyp];
                      shiftCoef   = shiftCoef + ntype[1 + ndtyp] * sizecoefs;  // shift coef   entre 2 types successif
                      shiftDonor = shiftDonor + ntype[1 + ndtyp];            // shift donor entre 2 types successif
                     }                                                        // type
                   count_rac += 1;
		  } // autorisation transfert
                }  // irac
            }      // pass_inst
  }    // omp

  delete [] RcvFields; delete [] DnrFields;

 #ifdef TimeShow
 #ifdef _OPENMP
    time_out = omp_get_wtime();
    ipt_timecount[2] = ipt_timecount[2] + time_out - time_in;
    time_in = omp_get_wtime();
 #endif
 #endif
    
    // send_buffer_time
#if defined(PROFILE_TRANSFERT)
#ifdef TimeShow
#ifdef _OPENMP
  double beg = omp_get_wtime();
#endif
#endif
#endif
  
#if defined(PROFILE_TRANSFERT)
#ifdef TimeShow
#ifdef _OPENMP
  double beg2 = omp_get_wtime();
#endif
#endif
#endif
  //if (has_data_to_send) { send_buffer.isend();  printf("envoi  ID  %d %d  %d \n", dest, TypeTransfert, nstep ); }
  if (has_data_to_send) { send_buffer.isend(); }
  //if (has_data_to_send) { send_buffer.isend(); send_buffer.test(); }

#ifdef TimeShow
#ifdef _OPENMP
    time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
#endif
#endif


#if defined(PROFILE_TRANSFERT)
#ifdef TimeShow
#ifdef _OPENMP
  double end = omp_get_wtime();
  isend_time += end - beg2;
  send_buffer_time += end - beg;
#endif
#endif
#endif
}

//=============================================================================
// Transfert de champs sous forme de numpy CMP lib GetData
// From zone
// Retourne une liste de numpy directement des champs interpoles
// in place + from zone + tc compact
//=============================================================================
void K_FASTC::getTransfersInter( E_Int& nbcom, E_Float**& ipt_ro, E_Int**& param_int, E_Float**& param_real, E_Int*& param_int_tc,
                                 std::pair<RecvQueue*, SendQueue*>*& pair_of_queue_loc, E_Float*& ipt_timecount)
 {
 
  if( nbcom != 0)
  {
     // Attente finalisation de la reception :
     assert(pair_of_queue_loc != NULL);

     vector<CMP::vector_view<E_Float> > recv_frp(4096);
     vector<E_Int> recv_nozone(4096);
     vector<E_Int> recv_nvarloc(4096);
     vector<E_Int> recv_size(4096);
     vector<CMP::vector_view<E_Int> > recv_listRc(4096);
     vector<CMP::vector_view<E_Int> > recv_fluxRc(4096);

     RecvQueue  pt_rcv_queue = *pair_of_queue_loc->first;



     while (not pt_rcv_queue.empty())
    {

#ifdef TimeShow
#ifdef _OPENMP
    E_Float time_in = omp_get_wtime();
#endif
#endif
       RecvQueue::iterator it = pt_rcv_queue.get_first_complete_message();
#ifdef TimeShow
#ifdef _OPENMP
    E_Float time_out = omp_get_wtime();
    ipt_timecount[6] = ipt_timecount[6] + time_out - time_in;
#endif
#endif

       if (it != pt_rcv_queue.end())
       {  // ok, une réception de prête*/

         CMP::RecvBuffer& recv_buffer = (*it).get_message_buffer();
#if defined(PROFILE_TRANSFERT)
#ifdef _OPENMP
         beg = omp_get_wtime();
#endif
#endif
         E_Int recv_nrac;
         recv_buffer >> recv_nrac;

         recv_frp.resize(recv_nrac);
         recv_nozone.resize(recv_nrac);
         recv_nvarloc.resize(recv_nrac);
         recv_listRc.resize(recv_nrac);
         recv_fluxRc.resize(recv_nrac);

         //recuperation des infos raccord en sequentiel
         for (E_Int irac = 0; irac < recv_nrac; ++irac)
         { 
          //recv_buffer >> recv_nozone[irac] >> recv_frp[irac] >> recv_listRc[irac];

          recv_buffer >> recv_nozone[irac] >> recv_nvarloc[irac] >> recv_frp[irac] >> recv_listRc[irac];
          //recv_buffer >> recv_nozone[irac] >> recv_nvarloc[irac] >> recv_frp[irac] >> recv_listRc[irac] >> recv_fluxRc[irac];

          E_Int nvars_loc = recv_nvarloc[irac]/100;
          E_Int nbFluCons = recv_nvarloc[irac] - nvars_loc*100;
          if(nbFluCons !=0){ recv_buffer >> recv_fluxRc[irac]; }
          recv_size[irac] = recv_listRc[irac].size();

          //recv_nvarloc[irac] = recv_frp[irac].size() / recv_size[irac];
           //printf("Nozone Recv Verif= %d , szrac: %d , irac: %d  \n", recv_nozone[irac],recv_size[irac], irac );
         }

#pragma omp parallel
         {
          for (E_Int irac = 0; irac < recv_nrac; ++irac)
          { 
            E_Int NoR = recv_nozone[irac];

            if (NoR == -999) continue;

            E_Int ilistrecv;
            E_Int nbRcvPts  = recv_size[irac];
            E_Int nvars_loc   = recv_nvarloc[irac]/100;
            E_Int nbFluCons = recv_nvarloc[irac] - nvars_loc*100;
            E_Int decal = param_int[ NoR ] [ NDIMDX ];
            if (NoR < 0) decal=0;

            E_Int shift_fluR = 0;
            for (E_Int nbflu = 0; nbflu < nbFluCons; nbflu++)
             {
             // stokage dans flux grossier pour nearmatch conservatif valable pour 5 et 6 Eq
              #include "FastC/Com/flux_conservatif_5eq_Recep.cpp"
             }


            //printf("Nozone Verif= %d nvar: %d , nflux: %d\n", NoR,nvars_loc, nbFluCons  );
            //fflush(stdout);
            if ( nvars_loc == 5)
            {
              #pragma omp for nowait
              for (E_Int irecv = 0; irecv <  nbRcvPts; ++irecv) 
               {

                ilistrecv = recv_listRc[irac] [irecv];

                //if (NoR==0 && irecv==0) printf("getTrans %f %d %d %d %d %p \n",recv_frp[irac][irecv], irecv, ilistrecv, irac , nbRcvPts, recv_frp[irac]);

                ipt_ro[NoR][ilistrecv          ] = recv_frp[irac][irecv];
                ipt_ro[NoR][ilistrecv +   decal] = recv_frp[irac][irecv + 1 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 2*decal] = recv_frp[irac][irecv + 2 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 3*decal] = recv_frp[irac][irecv + 3 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 4*decal] = recv_frp[irac][irecv + 4 * nbRcvPts];
               }
            } 
            else if ( nvars_loc == 6)
            {
            #pragma omp for nowait
            for (E_Int irecv = 0; irecv < nbRcvPts; ++irecv)
               {

                ilistrecv = recv_listRc[irac] [irecv]; 

                ipt_ro[NoR][ilistrecv          ] = recv_frp[irac][irecv];
                ipt_ro[NoR][ilistrecv +   decal] = recv_frp[irac][irecv + 1 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 2*decal] = recv_frp[irac][irecv + 2 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 3*decal] = recv_frp[irac][irecv + 3 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 4*decal] = recv_frp[irac][irecv + 4 * nbRcvPts];
                ipt_ro[NoR][ilistrecv + 5*decal] = recv_frp[irac][irecv + 5 * nbRcvPts];
               }
            } 
            else  {
                 for (int eq = 0; eq < nvars_loc; ++eq)
                 {
                  #pragma omp for nowait
                  for (E_Int irecv = 0; irecv < nbRcvPts; ++irecv)
                   {
                    ilistrecv = recv_listRc[irac] [irecv]; 

                    ipt_ro[NoR][ilistrecv + eq*decal] = recv_frp[irac][irecv + eq * nbRcvPts];
                   }  
                 }
                }
          }// rac
         }//  parallel....

#if defined(PROFILE_TRANSFERT)
#ifdef _OPENMP
         end = omp_get_wtime();
         recv_buffer_time += end - beg;
#endif
#endif
         pt_rcv_queue.pop(it);
       }  // end if (it != pt_msg_manager->end()infos
    }    // End  while (not pt_msg_manager->empty())
  }    // End if  nbcom=0 

  //SendQueue* pt_snd_queue =  pair_of_queue_loc->second;
  //pt_snd_queue->waitAll();
  //pt_snd_queue->clear();

#if defined(PROFILE_TRANSFERT)
  std::cout << "Temps actuel passe :  " << send_buffer_time << "( " << isend_time
            << " ) vs " << recv_buffer_time << std::endl;
#endif
}
#endif
