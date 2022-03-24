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
#include "FastS/fastS.h"
#include "FastS/param_solver.h"

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

#ifdef _MPI

#if defined(PROFILE_TRANSFERT)
static double send_buffer_time = 0.;
static double recv_buffer_time = 0.;
static double isend_time = 0.;
#endif

#ifdef _MPI
static std::pair<RecvQueue*, SendQueue*>* pair_of_queue = NULL;
static std::pair<RecvQueue*, SendQueue*>* pair_of_queue_IBC = NULL;
#endif

void K_FASTS::init_TransferInter(
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue) {
  int szCom;
  MPI_Comm_size(MPI_COMM_WORLD, &szCom);
  RecvQueue* pt_RecvQueue = new RecvQueue(szCom);
  SendQueue* pt_SendQueue = new SendQueue(szCom);
  pair_of_queue = new std::pair<RecvQueue*, SendQueue*>(pt_RecvQueue, pt_SendQueue);
//      std::make_pair(pt_RecvQueue, pt_SendQueue));
}

void K_FASTS::del_TransferInter(
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue) {
  delete pair_of_queue->first;
  delete pair_of_queue->second;
  delete pair_of_queue;
}
#endif

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
void K_FASTS::setInterpTransfersFastS(
  E_Float**& iptro_tmp, E_Int& vartype         , E_Int*& param_int_tc, E_Float*& param_real_tc , E_Int**& param_int     , E_Float**& param_real,
  E_Int*& linelets_int, E_Float*& linelets_real, E_Int& it_target    , E_Int& nidom            , E_Float*& ipt_timecount, E_Int& mpi,
  E_Int& nstep        , E_Int& nitmax          , E_Int& rk           , E_Int& exploc           , E_Int& numpassage)

{
  E_Int rank = 0;
  E_Int dest = 0;
  E_Int count= 0;


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
    E_Int TypeTransfert = 2;    
    E_Int nbcomIBC = param_int_tc[1];
    E_Int nbcomID  = param_int_tc[2+nbcomIBC];
    E_Int shift_graph = nbcomIBC + nbcomID + 2;

    #ifdef _MPI
    //std::pair<RecvQueue*, SendQueue*>* pair_of_queue;
    RecvQueue* pt_rcv_queue     = NULL;
    RecvQueue* pt_rcv_queue_IBC = NULL;


    E_Int nbcomIBC_nstep; E_Int pt_deb; E_Int nbcomID_nstep; 
    //if (rk ==3 and exploc==2 and mpi)
    if (exploc==1 and mpi)
      {
        pt_deb         = param_int_tc[1+nstep+(numpassage-1)*nitmax] + 2;
        nbcomIBC_nstep = param_int_tc[ pt_deb ];
	
      }
    else
      { pt_deb         = 1;
        nbcomIBC_nstep = nbcomIBC;
      }

    if (mpi)
    {
       
      if (pair_of_queue == NULL)
      {
         K_FASTS::init_TransferInter(pair_of_queue);
         K_FASTS::init_TransferInter(pair_of_queue_IBC);
      }

      #ifdef TimeShow
       time_in = omp_get_wtime();
      #endif

      pt_rcv_queue     = pair_of_queue->first;
      pt_rcv_queue_IBC = pair_of_queue_IBC->first;

      if (pt_rcv_queue_IBC->size() == 0)
        for (E_Int ircv = 1; ircv < nbcomIBC_nstep +1; ++ircv)
          {
            pt_rcv_queue_IBC->emplace_back(param_int_tc[pt_deb + ircv], 404);
            CMP::RecvBuffer& recv_buffer = pt_rcv_queue_IBC->back_message_buffer();
            recv_buffer.irecv();
          }
      else{
         assert(pt_rcv_queue_IBC->size() == nbcomIBC_nstep );
         for ( auto iterBuf = pt_rcv_queue_IBC->begin(); iterBuf != pt_rcv_queue_IBC->end(); ++iterBuf )
           {
            CMP::RecvBuffer& recv_buffer = iterBuf->get_message_buffer();
            recv_buffer.irecv();
           }
          }

      //cout << "proc : " << rank << "  nbcom : " << nbcomIBC_nstep << endl;

      #ifdef TimeShow
       E_Float time_out = omp_get_wtime();
       ipt_timecount[0] = ipt_timecount[0] + time_out -time_in;
      #endif

      int nb_send_buffer = 0;
      for (E_Int ip2p = 1; ip2p < param_int_tc[0] +1; ++ip2p)
      {
        E_Int ech  = param_int_tc[ip2p + shift_graph];
        dest       = param_int_tc[ech];

        if (dest != rank)  // Inter Process
        {
          nb_send_buffer += 1;
          TypeTransfert = 1;  
          #ifdef _MPI      
          K_FASTS::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
  	       			           param_int    , param_real   , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                           pair_of_queue_IBC, ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, nb_send_buffer); 
          #endif
        }
      }//loop comm p2p
    }  //mpi
    // _MPI
    #endif

    //comm local pour recouvrememnt
    for (E_Int ip2p = 1; ip2p < param_int_tc[0]+1; ++ip2p)
    {
      E_Int ech  = param_int_tc[ip2p+shift_graph];
      dest       = param_int_tc[ech];
      if (dest == rank)  // Intra Process
      { TypeTransfert = 1;
	K_FASTS::setInterpTransfersIntra(iptro_tmp    , vartype    , param_int_tc, param_real_tc,
                                         param_int    , param_real , linelets_int, linelets_real, TypeTransfert, it_target, nidom, ip2p, 
                                         ipt_timecount, nstep      , nitmax      , rk           , exploc       , numpassage); 
      }
    } //loop ip2p

    #ifdef TimeShow
      time_in = omp_get_wtime();
    #endif

    #ifdef _MPI
    if (mpi)
    {

      if (nbcomIBC_nstep != 0)

	{

	  //comm multi processus: wait + remplissage
	  K_FASTS::getTransfersInter(iptro_tmp, param_int, param_int_tc , pair_of_queue_IBC);

          #ifdef TimeShow
	  time_out = omp_get_wtime();
	  ipt_timecount[4] = ipt_timecount[4] + time_out -time_in;
	  time_in= omp_get_wtime();
          #endif

	}
      


      //if (rk ==3 and exploc==2 and mpi)
      if (exploc==1 and mpi)
      {
       pt_deb         = param_int_tc[2+nbcomIBC+nstep+(numpassage-1)*nitmax] + 2 +nbcomIBC+1;
       nbcomID_nstep  = param_int_tc[ pt_deb ];

      }
      else
      { pt_deb        = nbcomIBC + 2;
        nbcomID_nstep = nbcomID;
      }

      
     
      
      if (pt_rcv_queue->size() == 0 )
	for (E_Int ircv = 1; ircv < nbcomID_nstep +1; ++ircv)
	  {
	    pt_rcv_queue->emplace_back( param_int_tc[ pt_deb + ircv], 404);
	    CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
	    recv_buffer.irecv();
	  }     
      else {
        assert(pt_rcv_queue->size() == nbcomID_nstep );
        for ( auto iterBuf = pt_rcv_queue->begin(); iterBuf != pt_rcv_queue->end(); ++iterBuf )
        {
           CMP::RecvBuffer& recv_buffer = iterBuf->get_message_buffer();
           recv_buffer.irecv();
        }
       }

      #ifdef TimeShow
       time_out = omp_get_wtime();
       ipt_timecount[0] = ipt_timecount[0] + time_out -time_in;  
      #endif

      int cpt_send_buffer = 0;
      for (E_Int ip2p = 1; ip2p < param_int_tc[0]+1; ++ip2p)
        {
         E_Int ech  = param_int_tc[ip2p+shift_graph];
         dest       = param_int_tc[ech];

	 //cout << ech <<"  " << dest << endl;
    
         if (dest != rank)  // Inter Process
          {             
            TypeTransfert = 0;
            cpt_send_buffer += 1;
            K_FASTS::setInterpTransfersInter(iptro_tmp    , vartype      , param_int_tc, param_real_tc,
  	        	                     param_int    , param_real   , linelets_int, linelets_real, TypeTransfert, it_target , nidom , ip2p, 
                                             pair_of_queue, ipt_timecount, nstep       , nitmax       , rk           , exploc    , numpassage, cpt_send_buffer); 
          }
        }//loop ip2p
    } // Endif MPI
    #endif 

    //comm local pour recouvrememnt
    for (E_Int ip2p = 1; ip2p < param_int_tc[0]+1; ++ip2p)
    {
      E_Int ech  = param_int_tc[ip2p+shift_graph];
      dest       = param_int_tc[ech];
      if (dest == rank)  // Intra Process
      {        
        TypeTransfert = 0;
        K_FASTS::setInterpTransfersIntra(iptro_tmp, vartype   , param_int_tc, param_real_tc ,
                                         param_int, param_real, linelets_int, linelets_real, TypeTransfert , it_target, nidom, ip2p, ipt_timecount, nstep, nitmax, rk, exploc, numpassage); 
      }    
    }

    #ifdef TimeShow
     time_in = omp_get_wtime();
    #endif

    #ifdef _MPI
    if (mpi)
    {

      if (nbcomID_nstep != 0)

	{

       //comm multi processus: wait + remplissage
       K_FASTS::getTransfersInter(iptro_tmp, param_int, param_int_tc , pair_of_queue);



       #ifdef TimeShow
        time_out         = omp_get_wtime();
        ipt_timecount[4] = ipt_timecount[4] + time_out -time_in;

        outputfile << "Time in getTransfersInter "     << ipt_timecount[4] << std::endl;
        outputfile << "Time InterpTransfert (Intra)  " << ipt_timecount[1] << std::endl;
        outputfile << "Time in MPI send_buffer, irecv "<< ipt_timecount[0] << std::endl;
        outputfile << "Time InterpTransfert (Inter)  " << ipt_timecount[2] << std::endl;
        outputfile << "Nb com. p2p " << param_int_tc[0] +1 << std::endl;
        outputfile << std::endl << std::endl;
        outputfile.close();

        time_in = omp_get_wtime();
       #endif

	}
      
    } // MPI Second part (InterCOM ID)
    #endif      

  } //if  param_int_tc != Null
}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
void K_FASTS::setInterpTransfersIntra(

    E_Float**& ipt_ro, E_Int& varType, E_Int*& ipt_param_int_tc,
    E_Float*& ipt_param_real_tc, E_Int**& param_int, E_Float**& param_real, E_Int*& linelets_int, E_Float*& linelets_real,
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
  else
    nvars = 6;
  
  E_Int* ipt_cnd = NULL;  // ONLY FOR STRUCTURED
  
  E_Int nbcomIBC = ipt_param_int_tc[1];
  E_Int nbcomID  = ipt_param_int_tc[2 + nbcomIBC];
  
  E_Int shift_graph = nbcomIBC + nbcomID + 2;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech           = ipt_param_int_tc[NoTransfert + shift_graph];
  E_Int nrac          = ipt_param_int_tc[ech + 1];  // nb total de raccord
  E_Int nrac_inst     = ipt_param_int_tc[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel     = ipt_param_int_tc[ech + 3];  // nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady   = nrac - nrac_inst;        // nb total de raccord stationnaire

  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0; 
  E_Int pass_inst_fin=1;

  if (nrac_inst > 0) pass_inst_fin=2;

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(  size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation]; // Pour l explicite local

 
  // on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx = 0; E_Int ibcTypeMax=0;
  for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++) 
  {
    E_Int irac_deb = 0;
    E_Int irac_fin = nrac_steady;
    if(pass_inst == 1){ irac_deb = ipt_param_int_tc[ ech + 4 + it_target             ]; 
                        irac_fin = ipt_param_int_tc[ ech + 4 + it_target + timelevel ];
                      }

   for (E_Int irac = irac_deb; irac < irac_fin; irac++) 
    {
      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
      E_Int ibcType   = ipt_param_int_tc[shift_rac + nrac * 3];  

      if (ipt_param_int_tc[shift_rac + nrac * 10 + 1] > nbRcvPts_mx){ nbRcvPts_mx = ipt_param_int_tc[shift_rac + nrac * 10 + 1];}
      if( ibcType > ibcTypeMax)  ibcTypeMax =  ibcType;

      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      if(rk==3 and exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones en fonction de la ss-ite courante
	{
	  E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
	  E_Int levelD = ipt_param_int_tc[debut_rac + 25];
	  E_Int levelR = ipt_param_int_tc[debut_rac + 24];
	  E_Int cyclD  = nssiter/levelD;

	  //cout << levelD <<" " <<levelR << endl;

	  // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse   
	  if (levelD > levelR and num_passage == 1)		
	    {

	      if (nstep%cyclD==cyclD-1 or nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1)
		{
		  autorisation_transferts[pass_inst][irac_auto]=1;
		  if (ipt_param_int_tc[shift_rac + nrac * 10 + 1] > nbRcvPts_mx){ nbRcvPts_mx = ipt_param_int_tc[shift_rac + nrac * 10 + 1];}
		}

	    }
	  // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	  else if (levelD < levelR and num_passage == 1) 
	    {
	      if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)

		{
		  autorisation_transferts[pass_inst][irac_auto]=1;
		  if (ipt_param_int_tc[shift_rac + nrac * 10 + 1] > nbRcvPts_mx) {nbRcvPts_mx = ipt_param_int_tc[shift_rac + nrac * 10 + 1];}
		}
     

	    }
	  // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	  else if (levelD == levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1)

		{
		  autorisation_transferts[pass_inst][irac_auto]=1;
		  if (ipt_param_int_tc[shift_rac + nrac * 10 + 1] > nbRcvPts_mx){nbRcvPts_mx = ipt_param_int_tc[shift_rac + nrac * 10 + 1];}
		}
	      //else {continue;}

	    }
	  // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)   
	  else if (levelD == levelR and num_passage == 2)
	    {

	      //if (nstep%8 == 6){autorisation_transferts[pass_inst][irac]=1;}
	      if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1)
		{
		  autorisation_transferts[pass_inst][irac_auto]=1;
		  if (ipt_param_int_tc[shift_rac + nrac * 10 + 1] > nbRcvPts_mx){nbRcvPts_mx = ipt_param_int_tc[shift_rac + nrac * 10 + 1];}		
		}
	      //else {continue;}

	    }
	}
      else // Sinon, on autorise les transferts entre ttes les zones a ttes les ss-ite
	{
	  autorisation_transferts[pass_inst][irac_auto]=1;
	  if (ipt_param_int_tc[shift_rac + nrac * 10 + 1] > nbRcvPts_mx){nbRcvPts_mx = ipt_param_int_tc[shift_rac + nrac * 10 + 1];}
 
	}	
    }
  }

  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 14 * threadmax_sdm);
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
        // printf("ipass_inst = %d, level= %d \n",  ipass_inst, nrac_inst_level
        // );
        E_Int irac_deb = 0;
        E_Int irac_fin = nrac_steady;
        if (pass_inst == 1) {
          irac_deb = ipt_param_int_tc[ech + 4 + it_target];
          irac_fin = ipt_param_int_tc[ech + 4 + it_target + timelevel];
        }

        // printf("iracdeb=  %d, iracfin= %d \n", irac_deb, irac_fin  );
        for (E_Int irac = irac_deb; irac < irac_fin; irac++) {
	  
	 E_Int irac_auto= irac-irac_deb;
	 if (autorisation_transferts[pass_inst][irac_auto]==1)
	  {

          E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
          //printf("ipass_typ = %d, pass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,pass_inst,irac , ithread ); 
          // ipass_typ,ipass_inst,irac , ithread );
          E_Int ibcType = ipt_param_int_tc[shift_rac + nrac * 3];
          E_Int ibc = 1;
          if (ibcType < 0) ibc = 0;
          if (1 - ibc != ipass_typ) continue;

          E_Int NoD       = ipt_param_int_tc[shift_rac + nrac * 5     ];
          E_Int loc       = ipt_param_int_tc[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
          E_Int NoR       = ipt_param_int_tc[shift_rac + nrac * 11 + 1];
          E_Int nvars_loc = ipt_param_int_tc[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
          E_Int rotation  = ipt_param_int_tc[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal

          E_Float Pr    = param_real[ NoD ][ PRANDT ];
          E_Float Ts    = param_real[ NoD ][ TEMP0 ];
          E_Float Cs    = param_real[ NoD ][ CS ];
          E_Float muS   = param_real[ NoD ][ XMUL0 ];
          E_Float cv    = param_real[ NoD ][ CVINF ];
          E_Float gamma = param_real[ NoD ][ GAMMA ];

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

          E_Int nbRcvPts = ipt_param_int_tc[shift_rac + nrac * 10 + 1];
          // E_Int nbDonPts = ipt_param_int_tc[ shift_rac                ];

          E_Int pos;
          pos = ipt_param_int_tc[shift_rac + nrac * 7];      E_Int* ntype    = ipt_param_int_tc  + pos;
          pos = pos + 1 + ntype[0];                          E_Int* types    = ipt_param_int_tc  + pos;
          pos = ipt_param_int_tc[shift_rac + nrac * 6];      E_Int* donorPts = ipt_param_int_tc  + pos;
          pos = ipt_param_int_tc[shift_rac + nrac * 12 + 1]; E_Int* rcvPts   = ipt_param_int_tc  + pos;  // donor et receveur inverser car storage donor
          pos = ipt_param_int_tc[shift_rac + nrac * 8];    E_Float* ptrCoefs = ipt_param_real_tc + pos;

          E_Int nbInterpD = ipt_param_int_tc[shift_rac + nrac];
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
                        if (varType == 1 || varType == 11)
                          K_CONNECTOR::setIBCTransfersCommonVar1(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                                 xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                                                 xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                                                 xPI, xPI+nbRcvPts, xPI+nbRcvPts*2, 
                                                                 densPtr, densPtr+nbRcvPts, //dens + press
                                                                 densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                                                 densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                                                 densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,  
                                                                 ipt_tmp, size,
                                                                 gamma, cv, muS, Cs, Ts, Pr,
                                                                 vectOfDnrFields, vectOfRcvFields);
                        else if (varType == 2 || varType == 21)
                        {                                                                                      
                         K_CONNECTOR::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                                 xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                                                 xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                                                 xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                                                 densPtr, densPtr+nbRcvPts, //dens + press
                                                                 densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                                                 densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                                                 densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,  
                                                                 ipt_tmp, size,
                                                                 param_real[ NoD ],
                                                                 //gamma, cv, muS, Cs, Ts, Pr,
                                                                 vectOfDnrFields, vectOfRcvFields
                                                                ,nbptslinelets, linelets, indexlinelets);

                        }
                        else if (varType == 3 || varType == 31)
                          K_CONNECTOR::setIBCTransfersCommonVar3(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                                 xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                                                 xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                                                 xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                                                 densPtr, densPtr+nbRcvPts, //dens + press
                                                                 densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                                                 densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                                                 densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,  
                                                                 ipt_tmp, size,
                                                                 gamma, cv, muS, Cs, Ts, Pr,
                                                                 vectOfDnrFields, vectOfRcvFields);
                      }//ibc          
                      //*
                      //        } //chunk
                      //*/
                      ideb       = ideb + ifin;
                      shiftCoef  = shiftCoef   +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
                      shiftDonor= shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif

                   }// type 
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
void K_FASTS::setInterpTransfersInter(
    E_Float**& ipt_ro   , E_Int& varType   , E_Int*& ipt_param_int_tc, E_Float*& ipt_param_real_tc,
    E_Int**& param_int  , E_Float**& param_real, E_Int*& linelets_int    , E_Float*& linelets_real, 
    E_Int& TypeTransfert, E_Int& it_target, E_Int& nidom, E_Int& NoTransfert,
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue,
    E_Float*& ipt_timecount                          ,
    E_Int& nstep, E_Int& nssiter, E_Int& rk, E_Int& exploc, E_Int& num_passage, E_Int& nb_send_buffer)

{

#ifdef TimeShow
    E_Float time_in = omp_get_wtime();
#endif    

  E_Int pass_deb, pass_fin;
  if (TypeTransfert == 0) {
    pass_deb = 1;
    pass_fin = 2;
  }  // ID
  else if (TypeTransfert == 1) {
    pass_deb = 0;
    pass_fin = 1;
  }  // IBCD
  else {
    pass_deb = 0;
    pass_fin = 2;
  }  // ALL

  E_Int* ipt_cnd = NULL;  // ONLY FOR STRUCTURED

  E_Int imd, jmd, imdjmd, kmd, nvars, ndimdxD;
  E_Float* iptroD;

  if (varType <= 3 && varType >= 1)
    nvars = 5;
  else
    nvars = 6;

  E_Int nbcomIBC = ipt_param_int_tc[1];
  E_Int nbcomID  = ipt_param_int_tc[2 + nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 2;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech  = ipt_param_int_tc[NoTransfert + shift_graph];
  E_Int nrac = ipt_param_int_tc[ech + 1];  // nb total de raccord
  E_Int nrac_inst = ipt_param_int_tc[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel = ipt_param_int_tc[ech + 3];  // nb de pas de temps stocker pour
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
   nrac_inst_level = ipt_param_int_tc[ech + 4 + it_target + timelevel] - ipt_param_int_tc[ech + 4 + it_target] + 1; 
  } 
  
  // on dimension tableau travail pour IBC et pour transfert
  // E_Int nrac_inst_level = ipt_param_int_tc[ech + 4 + it_target + timelevel] -
  //                         ipt_param_int_tc[ech + 4 + it_target] + 1;

  std::vector<E_Float*> frp(nrac_steady + nrac_inst_level);

  int dest = ipt_param_int_tc[ipt_param_int_tc[NoTransfert + shift_graph]]; 

  int rank;
#ifdef _MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#endif
  SendQueue* pt_snd_queue = pair_of_queue->second;
  if (pt_snd_queue->size() < nb_send_buffer)
  {
      pt_snd_queue->emplace_back(dest, 404);
  }

  // Préparation du buffer d'envoi :
  CMP::SendBuffer& send_buffer = *(*pt_snd_queue)[nb_send_buffer-1].message_buffer;//back_message_buffer();
  
  // A partir d'ici pour allouer les tableaux à remplir
  std::vector<CMP::SendBuffer::PackedData*> pck_data;

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
  E_Int debut_rac;
  E_Int cycl;

  
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++) 
    {
        E_Int irac_deb = 0;
        E_Int irac_fin = nrac_steady;
        if ( pass_inst == 1 ) 
        {
            irac_deb = ipt_param_int_tc[ech + 4 + it_target];
            irac_fin = ipt_param_int_tc[ech + 4 + it_target + timelevel];
        }

        for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
        {
        
	  E_Int irac_auto= irac-irac_deb;
	  autorisation_transferts[pass_inst][irac_auto]=0;	  


	  if(rk==3 and exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones en fonction de la ss-ite courante
	    {
	      E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;     
	      E_Int levelD = ipt_param_int_tc[debut_rac + 25];
	      E_Int levelR = ipt_param_int_tc[debut_rac + 24];
	      E_Int cyclD  = nssiter/levelD;

	      // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse   
	      if (levelD > levelR and num_passage == 1)		
		{
		  if (nstep%cyclD==cyclD-1 or nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1){autorisation_transferts[pass_inst][irac_auto]=1;}
		  else {continue;}
		}
	      // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	      else if (levelD < levelR and num_passage == 1) 
		{
		  if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)
                     {autorisation_transferts[pass_inst][irac_auto]=1;}
		  else {continue;}
		}
	      // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	      else if (levelD == levelR and num_passage == 1)
		{
		  //cout << "coucou " << endl;
		  if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1) { autorisation_transferts[pass_inst][irac_auto]=1; }
		  else {continue;}
		}
	      // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)   
	      else if (levelD ==  levelR and num_passage == 2)
		{
		  //if (nstep%8 == 6){autorisation_transferts[pass_inst][irac]=1;}
		  if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1){autorisation_transferts[pass_inst][irac_auto]=1;}
		  else {continue;}
		}
	      else {continue;} 
	    }
	  else // Sinon, on autorise les transferts entre ttes les zones a ttes les ss-ite
	    {
	      autorisation_transferts[pass_inst][irac_auto]=1;
	    }	
	  
	  //cout << "autorisation_transferts[pass_inst][irac]= " << autorisation_transferts[pass_inst][irac] << endl;
	  if (autorisation_transferts[pass_inst][irac_auto]==1)
	   { 
	      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	      E_Int nbRcvPts = ipt_param_int_tc[shift_rac + nrac * 10 + 1];
    
	      E_Int ibcType = ipt_param_int_tc[shift_rac + nrac * 3];

              if (ibcType > ibcTypeMax){ ibcTypeMax= ibcType;}
	      E_Int ibc = 1;
	      if (ibcType < 0) ibc = 0;
	      if (TypeTransfert == 0 && ibc == 1) 
		{
		  continue;
		} 
	      else if (TypeTransfert == 1 && ibc == 0) 
		{
		  continue;
		}
	      if (nbRcvPts > nbRcvPts_mx) nbRcvPts_mx = nbRcvPts;
	      has_data_to_send |= (TypeTransfert == ibc);
	      count_rac += 1;

	     }// autorisation transfert	
        } // irac
    }  // pas inst 

  cout << "proc =  " << rank << "  count_rac= " << count_rac << endl;
  //cout << has_data_to_send << endl;

E_Float time_out;
if (has_data_to_send) {
#ifdef TimeShow
 time_out = omp_get_wtime();
 ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
 time_in  = omp_get_wtime();
#endif
    

    send_buffer << count_rac;

    count_rac = 0;
    
    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) {
          irac_deb = ipt_param_int_tc[ech + 4 + it_target];
          irac_fin = ipt_param_int_tc[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
      {
	
	E_Int irac_auto= irac-irac_deb;
	if (autorisation_transferts[pass_inst][irac_auto]==1)
	{ 
    
	    E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
      
	    E_Int ibcType = ipt_param_int_tc[shift_rac + nrac * 3];
	    E_Int ibc = 1;
	    if (ibcType < 0) ibc = 0;

	    if (TypeTransfert == 0 && ibc == 1) {
	      continue;
	    } else if (TypeTransfert == 1 && ibc == 0) {
	      continue;
	    }
	    E_Int nbRcvPts  = ipt_param_int_tc[shift_rac + nrac * 10 + 1];
	    E_Int nvars_loc = ipt_param_int_tc[shift_rac + nrac * 13 + 1];  // flag raccord rans/LES
	    E_Int Nozone = ipt_param_int_tc[shift_rac + nrac * 11 + 1];

	    send_buffer << Nozone;  

	    pck_data.push_back(&send_buffer.push_inplace_array(nvars_loc * nbRcvPts *
							       sizeof(E_Float)));    
	    E_Int PtlistDonor  = ipt_param_int_tc[shift_rac + nrac * 12 + 1];
	    E_Int* ipt_listRcv = ipt_param_int_tc + PtlistDonor;

	    send_buffer << CMP::vector_view<E_Int>(ipt_listRcv, nbRcvPts);
	    count_rac += 1;
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
          irac_deb = ipt_param_int_tc[ech + 4 + it_target];
          irac_fin = ipt_param_int_tc[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
      {
	//cout << autorisation_transferts[pass_inst][irac] << " " << pass_inst <<" "<< irac << endl;
	E_Int irac_auto= irac-irac_deb;
	if (autorisation_transferts[pass_inst][irac_auto]==1)
	 { 
	    E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	    E_Int ibcType   = ipt_param_int_tc[shift_rac + nrac * 3];
	    E_Int ibc = 1;
	    if (ibcType < 0) ibc = 0;

	    if (TypeTransfert == 0 && ibc == 1) {
	      continue;
	    } else if (TypeTransfert == 1 && ibc == 0) {
	      continue;
	    }

	    frp[count_rac] = pck_data[count_rac]->data<E_Float>();
	    count_rac += 1;
	  }// autorisation transfert
      } // irac
    } // pass_inst   
  } // if has data to send

//cout << "coucou" << endl;

#ifdef TimeShow
    time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
    time_in  = omp_get_wtime();
#endif    

  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 14 * threadmax_sdm);
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

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc, nocf, indCoef, noi, sizecoefs, Nbchunk, imd, jmd, imdjmd;

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
                    irac_deb = ipt_param_int_tc[ech + 4 + it_target];
                    irac_fin = ipt_param_int_tc[ech + 4 + it_target + timelevel];
                }

                for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
                {
		  //cout << autorisation_transferts[pass_inst][irac] << " " << pass_inst <<" "<< irac  << endl;
		  E_Int irac_auto= irac-irac_deb;
		  if (autorisation_transferts[pass_inst][irac_auto]==1)
		   {
 
                   E_Int shift_rac = ech + 4 + timelevel*2 + irac;
                   E_Int ibcType = ipt_param_int_tc[shift_rac + nrac * 3];
                   E_Int ibc = 1;
                   if (ibcType < 0) ibc = 0;
                   if (1 - ibc != ipass_typ) continue;
         
                   E_Int NoD       = ipt_param_int_tc[shift_rac + nrac * 5     ];
                   E_Int loc       = ipt_param_int_tc[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
                   E_Int nbRcvPts  = ipt_param_int_tc[shift_rac + nrac * 10 + 1];
                   E_Int NoR       = ipt_param_int_tc[shift_rac + nrac * 11 + 1];
                   E_Int nvars_loc = ipt_param_int_tc[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
                   E_Int rotation  = ipt_param_int_tc[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal
         
                   E_Float Pr    = param_real[ NoD ][ PRANDT ];
                   E_Float Ts    = param_real[ NoD ][ TEMP0 ];
                   E_Float Cs    = param_real[ NoD ][ CS ];
                   E_Float muS   = param_real[ NoD ][ XMUL0 ];
                   E_Float cv    = param_real[ NoD ][ CVINF ];
                   E_Float gamma = param_real[ NoD ][ GAMMA ];

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
                   pos               = ipt_param_int_tc[shift_rac + nrac * 7];
                   E_Int* ntype      = ipt_param_int_tc + pos;
                   pos               = pos + 1 + ntype[0];
                   E_Int* types      = ipt_param_int_tc + pos;
                   pos               = ipt_param_int_tc[shift_rac + nrac * 6];
                   E_Int* donorPts   = ipt_param_int_tc + pos;
                   pos               = ipt_param_int_tc[shift_rac + nrac * 8];
                   E_Float* ptrCoefs = ipt_param_real_tc + pos;
                   // pos               = ipt_param_int_tc[shift_rac + nrac * 12 + 1];
                   // E_Int* rcvPts     = ipt_param_int_tc +  pos;
         
                   E_Int nbInterpD = ipt_param_int_tc[shift_rac + nrac];
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
                        E_Int addrlinelets   = linelets_int[count_racIBC + 3 ];
                        E_Float* linelets    = linelets_real + addrlinelets;
                        E_Int* indexlinelets = linelets_int + linelets_int[1]+1 + 3;
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
                  E_Int NoR = ipt_param_int_tc[shift_rac + nrac * 11 + 1];
                  //if (ipt_param_int_tc[ech]==0) printf("No rac= %d , NoR= %d, NoD= %d, Ntype= %d, ptdeb= %d, ptfin= %d, NptD= %d, neq= %d, skip= %d, rank= %d, dest= %d,  thread= %d\n",
                  //irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin , 
                  //ipt_param_int_tc[ shift_rac + nrac*10+1  ], ipt_param_int_tc[ shift_rac + nrac*13+1  ], ipt_param_int_tc[ shift_rac + nrac*15+1  ], 
                  //rank, ipt_param_int_tc[ ech  ], ithread );
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
                    if ( varType == 1 || varType == 11 )
                      K_CONNECTOR::setIBCTransfersCommonVar1(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, 
                                                             xPC, xPC + nbRcvPts, xPC + nbRcvPts * 2, 
                                                             xPW, xPW + nbRcvPts, xPW + nbRcvPts * 2, 
                                                             xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
                                                             densPtr, densPtr+nbRcvPts, //dens + press
                                                             densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                                             densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                                             densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,  
                                                             ipt_tmp, size, gamma, cv, muS, Cs,
                                                             Ts,  Pr,vectOfDnrFields, vectOfRcvFields );
                    else if ( varType == 2 || varType == 21 )
                      K_CONNECTOR::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, 
                                                             xPC, xPC + nbRcvPts, xPC + nbRcvPts * 2, 
                                                             xPW, xPW + nbRcvPts, xPW + nbRcvPts * 2, 
                                                             xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
                                                             densPtr, densPtr+nbRcvPts, //dens + press
                                                             densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                                             densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                                             densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                                             ipt_tmp, size,
                                                             param_real[ NoD ],
                                                             //gamma, cv, muS, Cs, Ts, Pr,
                                                             vectOfDnrFields, vectOfRcvFields,
                                                             nbptslinelets, linelets, indexlinelets );
                    else if ( varType == 3 || varType == 31 )
                      K_CONNECTOR::setIBCTransfersCommonVar3(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, 
                                                             xPC, xPC + nbRcvPts, xPC + nbRcvPts * 2, 
                                                             xPW, xPW + nbRcvPts, xPW + nbRcvPts * 2, 
                                                             xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
                                                             densPtr, densPtr+nbRcvPts, //dens + press
                                                             densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                                             densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                                             densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,   
                                                             ipt_tmp, size, gamma, cv, muS, Cs,
                                                             Ts, Pr, vectOfDnrFields, vectOfRcvFields );
                  }  // ibc
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
    time_out = omp_get_wtime();
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

  

  if (has_data_to_send) {
    send_buffer.isend();
   }
  
#ifdef TimeShow
    time_out = omp_get_wtime();
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
void K_FASTS::getTransfersInter( E_Float**& ipt_roD, E_Int**& param_int, E_Int*& ipt_param_int_tc, std::pair<RecvQueue*, SendQueue*>*& pair_of_queue) {
  
  // Attente finalisation de la réception :
  assert(pair_of_queue != NULL);

  //cout << "ap_getTransfers" << endl;

  vector<CMP::vector_view<E_Float> > recv_frp(2048);
  vector<E_Int> recv_nozone(2048);
  vector<E_Int> recv_nvarloc(2048);
  vector<CMP::vector_view<E_Int> > recv_listRc(2048);

  RecvQueue pt_rcv_queue = *pair_of_queue->first;
  SendQueue* pt_snd_queue = pair_of_queue->second;

  while (not pt_rcv_queue.empty()) {

    RecvQueue::iterator it = pt_rcv_queue.get_first_complete_message();
    if (it != pt_rcv_queue.end()) {  // ok, une réception de prête*/

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
      size_t sz;

      for (E_Int irac = 0; irac < recv_nrac; ++irac) {
        recv_buffer >> recv_nozone[irac] >> recv_frp[irac] >> recv_listRc[irac];

        std::size_t sz = recv_listRc[irac].size();
        recv_nvarloc[irac] = recv_frp[irac].size() / sz;

#pragma omp parallel
        {
          E_Int ilistrecv;
          //E_Int decal = ipt_ndimdxD[recv_nozone[irac]];
          E_Int decal = param_int[ recv_nozone[irac] ] [ NDIMDX ];

          if (recv_nvarloc[irac] == 5) {
#pragma omp for
            for (int irecv = 0; irecv < sz; ++irecv) {

              ilistrecv = recv_listRc[irac] [irecv]; 
 
              ipt_roD[recv_nozone[irac]][ilistrecv          ] = recv_frp[irac][irecv];
              ipt_roD[recv_nozone[irac]][ilistrecv +   decal] = recv_frp[irac][irecv + 1 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 2*decal] = recv_frp[irac][irecv + 2 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 3*decal] = recv_frp[irac][irecv + 3 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 4*decal] = recv_frp[irac][irecv + 4 * sz];
            }  // end for (int irecv
          } else if (recv_nvarloc[irac] == 6) {
#pragma omp for
            for (int irecv = 0; irecv < sz; ++irecv) {

              ilistrecv = recv_listRc[irac] [irecv]; 

              ipt_roD[recv_nozone[irac]][ilistrecv          ] = recv_frp[irac][irecv];
              ipt_roD[recv_nozone[irac]][ilistrecv +   decal] = recv_frp[irac][irecv + 1 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 2*decal] = recv_frp[irac][irecv + 2 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 3*decal] = recv_frp[irac][irecv + 3 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 4*decal] = recv_frp[irac][irecv + 4 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 5*decal] = recv_frp[irac][irecv + 5 * sz];
            }  // end for (int irecv
          }

        }  // end omp parallel
      }  // end for ( int irac = ....

#if defined(PROFILE_TRANSFERT)
      end = omp_get_wtime();
      recv_buffer_time += end - beg;
#endif
      pt_rcv_queue.pop(it);
    }  // end if (it != pt_msg_manager->end()infos
  }    // End  while (not pt_msg_manager->empty())

  //cout << "wait"<< endl;

  pt_snd_queue->waitAll();
  pt_snd_queue->clear();
#if defined(PROFILE_TRANSFERT)
  std::cout << "Temps actuel passe :  " << send_buffer_time << "( " << isend_time
            << " ) vs " << recv_buffer_time << std::endl;
#endif
}
#endif
