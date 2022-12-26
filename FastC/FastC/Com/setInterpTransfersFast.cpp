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

#include "connector.h"
#include "FastC/fastc.h"
#include "FastC/param_solver.h"

#ifdef _MPI
#include <mpi.h>
#include "CMP/include/pending_message_container.h"
#include "CMP/include/recv_buffer.hpp"
#include "CMP/include/send_buffer.hpp"
#include "setInterpTransfersD.h"
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

static E_Int source_flag[2048];

#ifdef _MPI

#if defined(PROFILE_TRANSFERT)
static double send_buffer_time = 0.;
static double recv_buffer_time = 0.;
static double isend_time = 0.;
#endif

static std::pair<RecvQueue*, SendQueue*>* pair_of_queue = NULL;
static std::pair<RecvQueue*, SendQueue*>* pair_of_queue_IBC = NULL;

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
  E_Float**& iptro_tmp, E_Int& vartype         , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int     , E_Float**& param_real, E_Int*& ipt_omp,
  E_Int*& linelets_int, E_Float*& linelets_real, E_Int& it_target    , E_Int& nidom            , E_Float*& ipt_timecount, E_Int& mpi,
  E_Int& nstep        , E_Int& nitmax          , E_Int& rk           , E_Int& exploc           , E_Int& numpassage)

{
  E_Int rank = 0;
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

   E_Float time_in = omp_get_wtime();
  #endif

 //Swap (call to setInterpTransfer)
  if ( (param_int_tc != NULL) && (param_real_tc != NULL))
  {
    E_Int TypeTransfert ;    
    E_Int sizecomIBC = param_int_tc[2];
    E_Int sizecomID  = param_int_tc[3+sizecomIBC];
    E_Int shift_graph = sizecomIBC + sizecomID + 3;

    //printf("VERIF %d %d init= %d , nstep = %d , mpi= %d \n", sizecomID, sizecomIBC, param_int_tc[0], nstep, mpi);

    #ifdef _MPI
    //std::pair<RecvQueue*, SendQueue*>* pair_of_queue;
    RecvQueue* pt_rcv_queue     = NULL;
    RecvQueue* pt_rcv_queue_IBC = NULL;

    //premier passage dans transfert couche C depuis  mise a plat de tc
    if (param_int_tc[0]==0 and mpi)
    { E_Int szCom;
      MPI_Comm_size(MPI_COMM_WORLD, &szCom);
      for (E_Int proc = 0; proc < szCom; ++proc){source_flag[proc]=0;}

      E_Int ech        = param_int_tc[1 + shift_graph];
      timelevel_tc     = param_int_tc[ech + 3];
 
      if (pair_of_queue  != NULL)
        {
          K_FASTC::del_TransferInter(pair_of_queue);
          K_FASTC::del_TransferInter(pair_of_queue_IBC);
          //printf("ALLO %d %d init= %d , nstep = %d , mpi= %d \n", sizecomID, sizecomIBC, param_int_tc[0], nstep, mpi);
        }
      //flag transfer initilisé. Remise a zero dans miseAplat.
      param_int_tc[0]=1;
    }

    E_Int nbcomID_S; E_Int nbcomID_U;  E_Int nbcomIBC_S; E_Int pt_debID_S; E_Int pt_debID_U; E_Int pt_debIBC_S;

    // info Comm ID instationnaire et dtloc G Jeanmass
    E_Int iter = 1;
    if (exploc==1 and mpi) { iter = nstep; nbcomID_U =0;}
    else if(timelevel_tc !=0)
       {  
          pt_debID_U = param_int_tc[3+sizecomIBC+ iter + 1 + it_target] + 3 +sizecomIBC+1;
          nbcomID_U  = param_int_tc[ pt_debID_U ];
       }
    else { nbcomID_U =0;}


    // info Comm ID stationnaire
    pt_debID_S = param_int_tc[3+sizecomIBC+iter] + 3 +sizecomIBC+1;
    nbcomID_S  = param_int_tc[ pt_debID_S ];


    // info Comm IBC stationnaire
    pt_debIBC_S= param_int_tc[2+iter] + 3;
    nbcomIBC_S = param_int_tc[ pt_debIBC_S];

    if (mpi and (nbcomID_S != 0 or nbcomIBC_S != 0 or nbcomID_U != 0))
    {
      if (pair_of_queue     == NULL) { K_FASTC::init_TransferInter(pair_of_queue    );}
      if (pair_of_queue_IBC == NULL) { K_FASTC::init_TransferInter(pair_of_queue_IBC);}

      #ifdef TimeShow
       time_in = omp_get_wtime();
      #endif

    //
    //
    //Debut Transfert IBC
    //
    //
      pt_rcv_queue_IBC = pair_of_queue_IBC->first;

      if (pt_rcv_queue_IBC->size() == 0)
        {
          for (E_Int ircv = 1; ircv < nbcomIBC_S +1; ++ircv)
            {
              pt_rcv_queue_IBC->emplace_back(param_int_tc[pt_debIBC_S + ircv], 404);
              CMP::RecvBuffer& recv_buffer = pt_rcv_queue_IBC->back_message_buffer();
              recv_buffer.irecv();
              //printf("receptionIBM source  nstep= %d , tag= %d  size= %d \n",  nstep,  recv_buffer.tag(), recv_buffer.size());
            }
        }
      else
        {
          assert(pt_rcv_queue_IBC->size() == nbcomIBC_S );
          for ( auto iterBuf = pt_rcv_queue_IBC->begin(); iterBuf != pt_rcv_queue_IBC->end(); ++iterBuf )
            {
              CMP::RecvBuffer& recv_buffer = iterBuf->get_message_buffer();
              recv_buffer.irecv();
              
              //printf("receptionIBM OLD  nstep= %d , tag= %d  size= %d \n",  nstep, recv_buffer.tag(), recv_buffer.size());
            }
        }    

      #ifdef TimeShow
       E_Float time_out = omp_get_wtime();
       ipt_timecount[0] = ipt_timecount[0] + time_out -time_in;
      #endif

      //MPI_Barrier(MPI_COMM_WORLD);

      int nb_send_buffer = 0;
      for (E_Int ip2p = 1; ip2p < param_int_tc[1]+1; ++ip2p)
      {
        E_Int ech  = param_int_tc[ip2p + shift_graph];
        dest       = param_int_tc[ech];

        if (dest != rank)  // Inter Process ibc
        {
          nb_send_buffer += 1;
          TypeTransfert = 1;
          #ifdef _MPI
          //printf("inter IBM  %d %d \n", dest, nstep );
          K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
                                           param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                           pair_of_queue_IBC, ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer); 
          #endif
        }
      }//loop comm p2p
    }  //mpi
    // _MPI
    #endif

    //comm local pour recouvrememnt transfert IBC
    for (E_Int ip2p = 1; ip2p < param_int_tc[1]+1; ++ip2p)
    {
      E_Int ech  = param_int_tc[ip2p+shift_graph];
      dest       = param_int_tc[ech];
      if (dest == rank)  // Intra Process
      { TypeTransfert = 1;
        //printf("intra IBM  %d %d \n", dest, nstep );
        K_FASTC::setInterpTransfersIntra(iptro_tmp    , vartype    , param_int_tc, param_real_tc,
                                         param_int    , param_real , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target, nidom, ip2p, 
                                         ipt_timecount, nstep      , nitmax      , rk           , exploc       , numpassage); 
      }
    } //loop ip2p

    #ifdef TimeShow
      time_in = omp_get_wtime();
    #endif

    #ifdef _MPI
     if (mpi )
    {
      //comm multi processus: wait + remplissage IBC
        //printf("get IBM  %d %d %d \n", dest, nstep, nbcomIBC_S );
      //
      K_FASTC::getTransfersInter(nbcomIBC_S, iptro_tmp, param_int, param_int_tc , pair_of_queue_IBC);

      #ifdef TimeShow
       E_Float time_out = omp_get_wtime();
       ipt_timecount[4] = ipt_timecount[4] + time_out -time_in;
       time_in= omp_get_wtime();
      #endif
    }

    
    //MPI_Barrier(MPI_COMM_WORLD);
    //
    //Debut Transfert ID
    //
    //
    
    if (mpi and (nbcomID_S != 0 or nbcomID_U != 0) )
    {

/*
      pt_rcv_queue  = pair_of_queue->first;

      for (E_Int ircv = 1; ircv < nbcomID_S +1; ++ircv)
        {
         E_Int source =  param_int_tc[ pt_debID_S + ircv]
         if (source_flag[ source] == 0 )
          {
            pt_rcv_queue->emplace_back( source , 405);
            CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
            recv_buffer.irecv();

            source_flag[ source] = 1;
          }
         else
          {
          }
        }
*/
      pt_rcv_queue     = pair_of_queue->first;

      if (pt_rcv_queue->size() == 0 )
        {
          for (E_Int ircv = 1; ircv < nbcomID_S +1; ++ircv)
           {
            E_Int source = param_int_tc[ pt_debID_S + ircv];
            pt_rcv_queue->emplace_back( source , 405);
            CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
            recv_buffer.irecv();

            //printf("reception ID Steady source  %d %d \n", source, nstep );
           }
          for (E_Int ircv = 1; ircv < nbcomID_U +1; ++ircv)
           {
            E_Int source = param_int_tc[ pt_debID_U + ircv];
            pt_rcv_queue->emplace_back( source , 405);
            CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
            recv_buffer.irecv();

            //printf("reception ID Unsteady source  %d %d \n", source, nstep );
           }
        }
      else
        {   pt_rcv_queue->resize(nbcomID_S);
            assert(pt_rcv_queue->size() == nbcomID_S );
            for ( auto iterBuf = pt_rcv_queue->begin(); iterBuf != pt_rcv_queue->end(); ++iterBuf )
              {
                CMP::RecvBuffer& recv_buffer = iterBuf->get_message_buffer();
                recv_buffer.irecv();
                //printf("reception ID OLD     %d  \n", nstep );
              }
            for (E_Int ircv = 1; ircv < nbcomID_U +1; ++ircv)
              {
                 E_Int source = param_int_tc[ pt_debID_U + ircv];
                 pt_rcv_queue->emplace_back( source , 405);
                 CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
                 recv_buffer.irecv();

                //printf("reception ID Unsteady source  %d %d \n", source, nstep );
              }
        }


      #ifdef TimeShow
       E_Float time_out = omp_get_wtime();
       ipt_timecount[0] = ipt_timecount[0] + time_out -time_in;
      #endif

      int cpt_send_buffer = 0;
      for (E_Int ip2p = 1; ip2p < param_int_tc[1]+1; ++ip2p)
        {
         E_Int ech  = param_int_tc[ip2p+shift_graph];
         dest       = param_int_tc[ech];

          //printf("rank = %d, dest = %d \n",rank, dest );
         if (dest != rank)  // Inter Process Id
          {
            TypeTransfert = 0;
            cpt_send_buffer += 1;
            //printf(" inter ID  %d %d \n", dest, nstep );
            //
            K_FASTC::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
  	        	                     param_int    , param_real   , ipt_omp     , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                             pair_of_queue, ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, cpt_send_buffer); 
          }
        }//loop ip2p
    } // Endif MPI
    #endif


    //comm local pour recouvrememnt
    for (E_Int ip2p = 1; ip2p < param_int_tc[1]+1; ++ip2p)
    {
      E_Int ech  = param_int_tc[ip2p+shift_graph];
      dest       = param_int_tc[ech];
      if (dest == rank)  // Intra Process Id
      {
        TypeTransfert = 0;
        K_FASTC::setInterpTransfersIntra(iptro_tmp, vartype   , param_int_tc, param_real_tc ,
                                         param_int, param_real,  ipt_omp    ,linelets_int, linelets_real, TypeTransfert , it_target, nidom, ip2p, 
                                         ipt_timecount, nstep, nitmax, rk, exploc, numpassage); 
      }    
    }

    #ifdef TimeShow
     time_in = omp_get_wtime();
    #endif


    #ifdef _MPI
    if (mpi)
    {
       //comm multi processus: wait + remplissage ID
       //
       E_Int nbcomID = nbcomID_S + nbcomID_U;
       //printf(" get ID  %d %d \n",  nstep,  nbcomID);

       K_FASTC::getTransfersInter(nbcomID, iptro_tmp, param_int, param_int_tc , pair_of_queue);

       //printf(" apres get ID  %d %d \n", nstep,  nbcomID);

       #ifdef TimeShow
        E_Float time_out         = omp_get_wtime();
        ipt_timecount[4] = ipt_timecount[4] + time_out -time_in;

        outputfile << "Time in getTransfersInter "     << ipt_timecount[4] << std::endl;
        outputfile << "Time InterpTransfert (Intra)  " << ipt_timecount[1] << std::endl;
        outputfile << "Time in MPI send_buffer, irecv "<< ipt_timecount[0] << std::endl;
        outputfile << "Time InterpTransfert (Inter)  " << ipt_timecount[2] << std::endl;
        outputfile << "Nb com. p2p " << param_int_tc[1] +1 << std::endl;
        outputfile << std::endl << std::endl;
        outputfile.close();

        time_in = omp_get_wtime();
       #endif

    } // MPI Second part (InterCOM ID)
    #endif

    //K_FASTC::del_TransferInter(pair_of_queue);
    //K_FASTC::del_TransferInter(pair_of_queue_IBC);

  } //if  param_int_tc != Null


}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
void K_FASTC::setInterpTransfersIntra(

    E_Float**& ipt_ro, E_Int& varType, E_Int*& param_int_tc,
    E_Float*& param_real_tc, E_Int**& param_int, E_Float**& param_real,  E_Int*& ipt_omp, E_Int*& linelets_int, E_Float*& linelets_real,
    E_Int& TypeTransfert, E_Int& it_target, E_Int& nidom, E_Int& NoTransfert,
    E_Float*& ipt_timecount, E_Int& nstep, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& num_passage)
{

#ifdef TimeShow
  E_Float time_in = omp_get_wtime();
#endif

  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL

  // E_Int NoTransfert   = 1; // ONLY INTRA
  E_Int nvars;
  if (varType <= 3 && varType >= 1)
    nvars = 5;
  else if (varType == 4)
    nvars = param_int[0][NEQ_LBM];
  else
    nvars = 6;

  E_Int* ipt_cnd = NULL;  // ONLY FOR STRUCTURED
  
  E_Int sizecomIBC = param_int_tc[2];
  E_Int sizecomID  = param_int_tc[3+sizecomIBC];
  E_Int shift_graph = sizecomIBC + sizecomID + 3;

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
      if      (TypeTransfert == 0 && ibc == 1) { continue; } 
      else if (TypeTransfert == 1 && ibc == 0) { continue; }


      // Si on est en explicit local, on va autoriser les transferts entre certaines zones en fonction de la ss-ite courante
      if(exploc == 1)
	{
	  E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
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
	      E_Int nbRcvPts = param_int_tc[shift_rac + nrac * 10 + 1];
    
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

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    // 1ere pass_typ: IBC
    // 2eme pass_typ: transfert
    //

    E_Int count_racIBC = 0;

     for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++)
     {
      // 1ere pass_inst: les raccord fixe
      // 2eme pass_inst: les raccord instationnaire
      for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
      {
        //  printf("pass %d %d %d %d \n", rank ,ipass_typ, pass_inst, it_target ); 
        // printf("ipass_inst = %d, level= %d \n",  ipass_inst, nrac_inst_level
        // );
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

          E_Int NoD       = param_int_tc[shift_rac + nrac * 5     ];
          E_Int NoR       = param_int_tc[shift_rac + nrac * 11 + 1];
          E_Int nvars_loc = param_int_tc[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
          E_Int rotation  = param_int_tc[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal

          //E_Float Pr    = param_real[ NoD ][ PRANDT ];
          //E_Float Ts    = param_real[ NoD ][ TEMP0 ];
          //E_Float Cs    = param_real[ NoD ][ CS ];
          //E_Float muS   = param_real[ NoD ][ XMUL0 ];
          //E_Float cv    = param_real[ NoD ][ CVINF ];
          //E_Float gamma = param_real[ NoD ][ GAMMA ];

          E_Int meshtype = param_int[ NoD ][ MESHTYPE ] ;

          E_Int cnNfldD = 0;
          E_Int* ptrcnd = NULL;


          

          for (E_Int eq = 0; eq < nvars_loc; eq++) {
            vectOfRcvFields[eq] = ipt_ro[NoR] + eq * param_int[NoR][ NDIMDX];
            vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX];
          }
          imd = param_int[ NoD ][NIJK  ];
          jmd = param_int[ NoD ][NIJK+1];
          // }

          imdjmd = imd * jmd;

          ////
          //  Interpolation parallele
          ////
          ////

          E_Int nbRcvPts = param_int_tc[shift_rac + nrac * 10 + 1];
          // E_Int nbDonPts = param_int_tc[ shift_rac                ];

          E_Int pos;
          pos = param_int_tc[shift_rac + nrac * 7];      E_Int* ntype    = param_int_tc  + pos;
          pos = pos + 1 + ntype[0];                      E_Int* types    = param_int_tc  + pos;
          pos = param_int_tc[shift_rac + nrac * 6];      E_Int* donorPts = param_int_tc  + pos;
          pos = param_int_tc[shift_rac + nrac * 12 + 1]; E_Int* rcvPts   = param_int_tc  + pos;  // donor et receveur inverser car storage donor
          pos = param_int_tc[shift_rac + nrac * 8];    E_Float* ptrCoefs = param_real_tc + pos;

          E_Int nbInterpD = param_int_tc[shift_rac + nrac];
          E_Float* xPC = NULL;
          E_Float* xPI = NULL;
          E_Float* xPW = NULL;
          E_Float* densPtr = NULL;
          E_Float* linelets    = NULL;
          E_Int* indexlinelets = NULL;
          E_Int nbptslinelets  = 0;
          if (ibc == 1) {
            xPC = ptrCoefs + nbInterpD;
            xPI = ptrCoefs + nbInterpD + 3 * nbRcvPts;
            xPW = ptrCoefs + nbInterpD + 6 * nbRcvPts;
            densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;

            if (linelets_int != NULL )
               {
               nbptslinelets        = linelets_int[0];
               E_Int addrlinelets   = linelets_int[count_racIBC + 3 ];
               linelets             = linelets_real + addrlinelets;
               indexlinelets        = linelets_int + linelets_int[1]+1 + 3;
               count_racIBC         = count_racIBC + 1;
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
            //         E_Int size_min=   16;
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

            // printf(" irac= %d, NoR= %d, nvar=  %d, NoD= %d, Rans=  %d, rot=
            // %d, fin= %d, ithread= %d \n", irac, NoR, nvars_loc, NoD,
            // ipass_inst ,rotation, pt_fin , ithread );  if(ithread <=8 &&
            // NoD==83 )  printf(" shift %d  %d %d %d  %d %d %d  %d \n", irac,
            // NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );
            // if(ithread <=8 && NoR==114 )  printf(" new   %d  %d %d %d  %d %d
            // %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  ,
            // type, ithread );

                      noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
                      indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
                      if     (nvars_loc==5)
                      {
            #           include "commonInterpTransfers_reorder_5eq.h"
                      }
                      else if(nvars_loc==6)
                      {
            #           include "commonInterpTransfers_reorder_6eq.h"
                      }
                      else
                      {
			//            #           include "LBM/commonInterpTransfers_reorder_neq.h"
                      }

                      // Prise en compte de la periodicite par rotation
                      if (rotation == 1)
                      {
                       E_Float* angle = ptrCoefs + nbInterpD;
            #          include "includeTransfers_rotation.h"
                      }
                      // ibc
                      if (ibc == 1)
			{
			  K_CONNECTOR::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                                 xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                                                 xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                                                 xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                                                 densPtr, 
                                                                 ipt_tmp, size,
                                                                 param_real[ NoD ],
                                                                 //gamma, cv, muS, Cs, Ts, Pr,
                                                                 vectOfDnrFields, vectOfRcvFields,
                                                                 nbptslinelets, linelets, indexlinelets);
                        
                      }//ibc
                      //*
                      //        } //chunk
                      //*/
                      ideb       = ideb + ifin;
                      shiftCoef  = shiftCoef   +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
                      shiftDonor= shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif

                   }// type
              #pragma omp barrier 
	          } //autorisation transfert
                }//irac
               }//pass_inst
              #pragma omp barrier
    }  // ipass
  }    // omp

#ifdef TimeShow
  E_Float time_out = omp_get_wtime();
  ipt_timecount[1] = ipt_timecount[1] + time_out - time_in;
  time_in = omp_get_wtime();
#endif

  delete[] ipt_cnd;
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
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue_loc,
    E_Float*& ipt_timecount                          ,
    E_Int& nstep, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& num_passage, E_Int& nb_send_buffer)

{
#ifdef TimeShow
    E_Float time_in = omp_get_wtime();
#endif

  E_Int pass_deb, pass_fin, etiquette;
  if (TypeTransfert == 0) {
    pass_deb = 1;
    pass_fin = 2;
    etiquette=405;
  }  // ID
  else if (TypeTransfert == 1) {
    pass_deb = 0;
    pass_fin = 1;
    etiquette=404;
  }  // IBCD
  else {
    pass_deb = 0;
    pass_fin = 2;
    etiquette=404;
  }  // ALL

  E_Int nvars;

  if (varType <= 3 && varType >= 1)
    nvars = 5;
  else if (varType == 4)
    nvars = param_int[0][NEQ_LBM];
  else
    nvars = 6;

  E_Int sizecomIBC = param_int_tc[2];
  E_Int sizecomID  = param_int_tc[3+sizecomIBC];
  E_Int shift_graph = sizecomIBC + sizecomID + 3;

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

#ifdef TimeShow
    time_in = omp_get_wtime();
#endif

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

  // A partir d'ici pour allouer les tableaux a� remplir
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
          E_Int ibcType = param_int_tc[shift_rac + nrac * 3];

          if (ibcType > ibcTypeMax){ ibcTypeMax= ibcType;}
          E_Int ibc = 1;
	  if (ibcType < 0) ibc = 0;
	  if      (TypeTransfert == 0 && ibc == 1) { continue; } 
	  else if (TypeTransfert == 1 && ibc == 0) { continue; }

	  if(exploc == 1)
	    {
	      E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
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
	      E_Int nbRcvPts = param_int_tc[shift_rac + nrac * 10 + 1];

	      if (nbRcvPts > nbRcvPts_mx) nbRcvPts_mx = nbRcvPts;
	      has_data_to_send |= (TypeTransfert == ibc);
	      //has_data_to_send = true;
	      count_rac += 1;

	     }// autorisation transfert
        } // irac
    }  // pas inst


if (has_data_to_send) {
#ifdef TimeShow
 E_Float time_out = omp_get_wtime();
 ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
 time_in  = omp_get_wtime();
#endif


    send_buffer << count_rac;

    //count_rac = 0;

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

	    E_Int nbRcvPts  = param_int_tc[shift_rac + nrac * 10 + 1];
	    E_Int nvars_loc = param_int_tc[shift_rac + nrac * 13 + 1];  // flag raccord rans/LES
	    E_Int Nozone    = param_int_tc[shift_rac + nrac * 11 + 1];

            // on determine un No zone pipeau pour skipper remplissage inutile en implicit local
            E_Int NoD       = param_int_tc[shift_rac + nrac * 5     ];
            E_Int Nozone_loc   = Nozone; 
            E_Int nbRcvPts_loc = nbRcvPts;
            if (impli_local[NoD] == 0) {Nozone_loc = -999; nbRcvPts_loc=1;}

            //if(Nozone==5 and count_rac==0 ) { printf("Nozone_loc %d \n", Nozone_loc);}

	    send_buffer << Nozone_loc;  

	    pck_data.push_back(&send_buffer.push_inplace_array(nvars_loc * nbRcvPts_loc * sizeof(E_Float) ));    
	    E_Int PtlistDonor  = param_int_tc[shift_rac + nrac * 12 + 1];
	    E_Int* ipt_listRcv = param_int_tc + PtlistDonor;

	    send_buffer << CMP::vector_view<E_Int>(ipt_listRcv, nbRcvPts_loc);

            //printf("size rac %d %d %d %d \n", nbRcvPts_loc*nvars_loc, count_rac, irac, pass_inst);

	    //count_rac += 1;
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
	    E_Int shift_rac = ech + 4 + timelevel * 2 + irac;

	    frp[count_rac] = pck_data[count_rac]->data<E_Float>();
	    count_rac += 1;
	  }// autorisation transfert
      } // irac
    } // pass_inst
  } // if has data to send

#ifdef TimeShow
    E_Float time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
    time_in  = omp_get_wtime();
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

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    // 1ere pass: IBC
    // 2eme pass: transfert
    //

    E_Int count_racIBC = 0;

    for ( E_Int ipass_typ = pass_deb; ipass_typ < pass_fin; ipass_typ++ )
        {
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
                  E_Int NoD       = param_int_tc[shift_rac + nrac * 5 ];
                  E_Int nvars_loc = param_int_tc[shift_rac + nrac * 13 +1];  // neq fonction raccord rans/LES

		  E_Int irac_auto = irac-irac_deb;

                  // on determine un envoi pipeau de taille 1*neq pour skipper remplissage inutile en implicit local
		  if ( impli_local[NoD]==0 and autorisation_transferts[pass_inst][irac_auto]==1)
                    {
                      //for (E_Int eq = 0; eq < nvars_loc; eq++) { frp[count_rac][eq] =0; }
                      count_rac += 1;
                    }

                  if (autorisation_transferts[pass_inst][irac_auto]==1 and impli_local[NoD]==1)
		   {

                     E_Int ibcType = param_int_tc[shift_rac + nrac * 3];

                     E_Int ibc = 1;
	             if (ibcType < 0) ibc = 0;

                     //E_Int loc       = param_int_tc[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
                     E_Int nbRcvPts  = param_int_tc[shift_rac + nrac * 10 + 1];
                     E_Int NoR       = param_int_tc[shift_rac + nrac * 11 + 1];
                     E_Int rotation  = param_int_tc[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal

                     //E_Float Pr    = param_real[ NoD ][ PRANDT ];
                     //E_Float Ts    = param_real[ NoD ][ TEMP0 ];
                     //E_Float Cs    = param_real[ NoD ][ CS ];
                     //E_Float muS   = param_real[ NoD ][ XMUL0 ];
                     //E_Float cv    = param_real[ NoD ][ CVINF ];
                     //E_Float gamma = param_real[ NoD ][ GAMMA ];

                     E_Int meshtype = 1;  // ONLY FOR STRUCTURE ipt_ndimdxD[NoD + nidom*6];
                     E_Int cnNfldD  = 0;
                     E_Int* ptrcnd  = NULL;

                     // printf("navr_loc %d %d %d \n", nvars_loc, nvars, Rans);
                     for (E_Int eq = 0; eq < nvars_loc; eq++) {
                       vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                       vectOfDnrFields[eq] = ipt_ro[NoD]    + eq * param_int[ NoD ][ NDIMDX ];
                     }

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
                     // pos               = param_int_tc[shift_rac + nrac * 12 + 1];
                     // E_Int* rcvPts     = param_int_tc +  pos;
  
                     E_Int nbInterpD = param_int_tc[shift_rac + nrac];
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
                      //E_Int NoR = param_int_tc[shift_rac + nrac * 11 + 1];
                      //if (param_int_tc[ech]==0) printf("No rac= %d , NoR= %d, NoD= %d, Ntype= %d, ptdeb= %d, ptfin= %d, NptD= %d, neq= %d, skip= %d, rank= %d, dest= %d,  thread= %d\n",
                      //irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin , 
                      //param_int_tc[ shift_rac + nrac*10+1  ], param_int_tc[ shift_rac + nrac*13+1  ], param_int_tc[ shift_rac + nrac*15+1  ], 
                      //rank, param_int_tc[ ech  ], ithread );
                      if ( nvars_loc == 5 ) {
#include "commonInterpTransfersD_reorder_5eq.h"
                        } else if ( nvars_loc == 6 ) {
#include "commonInterpTransfersD_reorder_6eq.h"
                        } else {
		    //#include "LBM/commonInterpTransfersD_reorder_neq.h"
                        }
                      // Prise en compte de la periodicite par rotation
                      if ( rotation == 1 ) {
                          E_Float* angle = ptrCoefs + nbInterpD;
#include "includeTransfersD_rotation.h"
                          }

                      // ibc
                      if (ibc == 1)
                      {
                        // tableau temporaire pour utiliser la routine commune K_CONNECTOR::setIBCTransfersCommon
                        for ( E_Int noind = pt_deb; noind < pt_fin; noind++ ) rcvPts[noind] = noind;
	  	        K_CONNECTOR::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, 
			                                      xPC, xPC + nbRcvPts, xPC + nbRcvPts * 2, 
							      xPW, xPW + nbRcvPts, xPW + nbRcvPts * 2, 
							      xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
							      densPtr, 
							      ipt_tmp, size,
							      param_real[ NoD ],
							      //gamma, cv, muS, Cs, Ts, Pr,
							      vectOfDnrFields, vectOfRcvFields,
							      nbptslinelets, linelets, indexlinelets);
                  
                      }  // ibc
	              E_Int PtlistDonor  = param_int_tc[shift_rac + nrac * 12 + 1];
	              E_Int* ipt_listRcv = param_int_tc + PtlistDonor;

                      //        } //chunk
                      ideb        = ideb + ifin;
                      shiftCoef   = shiftCoef + ntype[1 + ndtyp] * sizecoefs;  // shift coef   entre 2 types successif
                      shiftDonor = shiftDonor + ntype[1 + ndtyp];            // shift donor entre 2 types successif
                     }                                                        // type
                   count_rac += 1;
		  } // autorisation transfert
                }  // irac
            }      // pass_inst
#pragma omp barrier
        }  // ipass
  }    // omp

 #ifdef TimeShow
    E_Float time_out = omp_get_wtime();
    ipt_timecount[2] = ipt_timecount[2] + time_out - time_in;
    time_in = omp_get_wtime();
 #endif

    // send_buffer_time
#if defined(PROFILE_TRANSFERT)
  double beg = omp_get_wtime();
#endif

#if defined(PROFILE_TRANSFERT)
  double beg2 = omp_get_wtime();
#endif

  //if (has_data_to_send) { send_buffer.isend();  printf("envoi  ID  %d %d  %d \n", dest, TypeTransfert, nstep ); }
  if (has_data_to_send) { send_buffer.isend(); }

#ifdef TimeShow
    E_Float time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
#endif



#if defined(PROFILE_TRANSFERT)
  double end = omp_get_wtime();
  isend_time += end - beg2;
  send_buffer_time += end - beg;
#endif

}

//=============================================================================
// Transfert de champs sous forme de numpy CMP lib GetData
// From zone
// Retourne une liste de numpy directement des champs interpoles
// in place + from zone + tc compact
//=============================================================================
void K_FASTC::getTransfersInter( E_Int& nbcom, E_Float**& ipt_ro, E_Int**& param_int, E_Int*& param_int_tc, std::pair<RecvQueue*, SendQueue*>*& pair_of_queue_loc) {
 
  if( nbcom != 0)
  {
     // Attente finalisation de la reception :
     assert(pair_of_queue_loc != NULL);

     vector<CMP::vector_view<E_Float> > recv_frp(2048);
     vector<E_Int> recv_nozone(2048);
     vector<E_Int> recv_nvarloc(2048);
     vector<E_Int> recv_size(2048);
     vector<CMP::vector_view<E_Int> > recv_listRc(2048);

     RecvQueue  pt_rcv_queue = *pair_of_queue_loc->first;

     while (not pt_rcv_queue.empty())
    {

       RecvQueue::iterator it = pt_rcv_queue.get_first_complete_message();
       if (it != pt_rcv_queue.end())
       {  // ok, une réception de prête*/

         CMP::RecvBuffer& recv_buffer = (*it).get_message_buffer();
#if defined(PROFILE_TRANSFERT)
         beg = omp_get_wtime();
#endif
         E_Int recv_nrac;
         recv_buffer >> recv_nrac;

         recv_frp.resize(recv_nrac);
         recv_nozone.resize(recv_nrac);
         recv_nvarloc.resize(recv_nrac);
         recv_listRc.resize(recv_nrac);

         //recuperation des infos raccord en sequentiel
         for (E_Int irac = 0; irac < recv_nrac; ++irac)
         { 
          recv_buffer >> recv_nozone[irac] >> recv_frp[irac] >> recv_listRc[irac];

          recv_size[irac] = recv_listRc[irac].size();
          recv_nvarloc[irac] = recv_frp[irac].size() / recv_size[irac];
          
           //printf("Nozone Verif= %d %d %d  \n", recv_nozone[irac],recv_size[irac], irac );
         }

#pragma omp parallel
         {
          for (E_Int irac = 0; irac < recv_nrac; ++irac)
          { 
            E_Int NoR = recv_nozone[irac];
            //printf("Nozone Verif= %d \n", NoR);
            //fflush(stdout);

            if (NoR == -999) continue;

            E_Int ilistrecv;
            E_Int sz = recv_size[irac];
            E_Int decal = param_int[ NoR ] [ NDIMDX ];
            if (NoR < 0) decal=0;

            if (recv_nvarloc[irac] == 5)
            {
              //#pragma omp for nowait
              #pragma omp for
              for (size_t irecv = 0; irecv < sz; ++irecv) 
               {

                ilistrecv = recv_listRc[irac] [irecv];

                ipt_ro[NoR][ilistrecv          ] = recv_frp[irac][irecv];
                ipt_ro[NoR][ilistrecv +   decal] = recv_frp[irac][irecv + 1 * sz];
                ipt_ro[NoR][ilistrecv + 2*decal] = recv_frp[irac][irecv + 2 * sz];
                ipt_ro[NoR][ilistrecv + 3*decal] = recv_frp[irac][irecv + 3 * sz];
                ipt_ro[NoR][ilistrecv + 4*decal] = recv_frp[irac][irecv + 4 * sz];
               }
            } 
            else if (recv_nvarloc[irac] == 6)
            {
            //#pragma omp for nowait
            #pragma omp for
            for (E_Int irecv = 0; irecv < sz; ++irecv)
               {

                ilistrecv = recv_listRc[irac] [irecv]; 

                ipt_ro[NoR][ilistrecv          ] = recv_frp[irac][irecv];
                ipt_ro[NoR][ilistrecv +   decal] = recv_frp[irac][irecv + 1 * sz];
                ipt_ro[NoR][ilistrecv + 2*decal] = recv_frp[irac][irecv + 2 * sz];
                ipt_ro[NoR][ilistrecv + 3*decal] = recv_frp[irac][irecv + 3 * sz];
                ipt_ro[NoR][ilistrecv + 4*decal] = recv_frp[irac][irecv + 4 * sz];
                ipt_ro[NoR][ilistrecv + 5*decal] = recv_frp[irac][irecv + 5 * sz];
               }
            } 
            else  {
                 for (int eq = 0; eq < recv_nvarloc[irac]; ++eq)
                 {
                  #pragma omp for nowait
                  for (E_Int irecv = 0; irecv < sz; ++irecv)
                   {
                    ilistrecv = recv_listRc[irac] [irecv]; 

                    ipt_ro[NoR][ilistrecv + eq*decal] = recv_frp[irac][irecv + eq * sz];
                   }  
                 }
                }
          }// rac
         }//  parallel....

#if defined(PROFILE_TRANSFERT)
         end = omp_get_wtime();
         recv_buffer_time += end - beg;
#endif
         pt_rcv_queue.pop(it);
       }  // end if (it != pt_msg_manager->end()infos
    }    // End  while (not pt_msg_manager->empty())
  }    // End if  nbcom=0 

  SendQueue* pt_snd_queue =  pair_of_queue_loc->second;
  pt_snd_queue->waitAll();
  pt_snd_queue->clear();

#if defined(PROFILE_TRANSFERT)
  std::cout << "Temps actuel passe :  " << send_buffer_time << "( " << isend_time
            << " ) vs " << recv_buffer_time << std::endl;
#endif
}
#endif
