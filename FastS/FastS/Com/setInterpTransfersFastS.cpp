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
#include "fastS.h"

#ifdef _MPI
#include <mpi.h>
#include "CMP/include/pending_message_container"
#include "CMP/include/recv_buffer.hpp"
#include "CMP/include/send_buffer.hpp"
#include "setInterpTransfersD.h"
#endif
#include <utility>

using namespace std;
using namespace K_FLD;

E_Float time_in;
E_Float time_out;
E_Int timeShowFastS = 0;

#ifdef _MPI

#if defined(PROFILE_TRANSFERT)
static double send_buffer_time = 0.;
static double recv_buffer_time = 0.;
static double isend_time = 0.;
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
  E_Float**& iptro_tmp, E_Int*& ipt_ndimdx_trans, E_Int*& param_int_tc, E_Float*& param_real_tc ,
  E_Int*& param_bci, E_Float*& param_bcf, E_Int& it_target, E_Int& nidom, E_Float*& ipt_timecount)

{
  E_Int rank = 0;
  E_Int dest = 0;
  int init   = 0;
#ifdef _MPI
  // Check if MPI_Init has been called by converter...
  MPI_Initialized( &init ); 
  // ( init == 1 ) => it is an MPI run
  if(init)
  {
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);  
  }
#endif

  if (rank == 0 and timeShowFastS) {
    time_in = omp_get_wtime();
  }

 //Swap (call to setInterpTransfer)
  if ( (param_int_tc != NULL) && (param_real_tc != NULL))
  {
    E_Int TypeTransfert = 2;    
    E_Int nbcomIBC = param_int_tc[1];
    E_Int nbcomID  = param_int_tc[2+nbcomIBC];
    E_Int shift_graph = nbcomIBC + nbcomID + 2;

#ifdef _MPI

std::pair<RecvQueue*, SendQueue*>* pair_of_queue;

if (init)
{
#ifdef TimeShow
  if ( rank == 0)
  {
    time_init = omp_get_wtime();
  }
#endif 

  K_FASTS::init_TransferInter(pair_of_queue);

#ifdef TimeShow
  if(rank == 0)
  {
   time_init = omp_get_wtime();
  }
#endif

  RecvQueue* pt_rcv_queue = pair_of_queue->first;
  
  for (E_Int ircv = 1; ircv < nbcomIBC +1; ++ircv)
  {
   pt_rcv_queue->emplace_back(param_int_tc[ircv+1], 404);
   CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
   recv_buffer.irecv();
  }

#ifdef TimeShow
  if(rank == 0)
  {
    time_COM = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_COM -time_init;
  }
#endif

  for (E_Int ip2p = 1; ip2p < param_int_tc[0] +1; ++ip2p)
  {
    E_Int ech  = param_int_tc[ip2p + shift_graph];
    dest       = param_int_tc[ech];
    if (dest != rank)  // Inter Process
    {
      TypeTransfert = 1;        
      K_FASTS::setInterpTransfersInter(iptro_tmp,     ipt_ndimdx_trans,  param_int_tc,  param_real_tc ,
                                           param_bci, param_bcf, TypeTransfert , it_target, nidom, ip2p, pair_of_queue, ipt_timecount); 
    }
    else
    {
		TypeTransfert = 1;
    K_FASTS::setInterpTransfersIntra(iptro_tmp,     ipt_ndimdx_trans,  param_int_tc,  param_real_tc ,
                                             param_bci, param_bcf, TypeTransfert , it_target, nidom, ip2p, ipt_timecount); 
    }
    
  }

#ifdef TimeShow
  if( rank == 0 )
  {
    time_init = omp_get_wtime();
  }
#endif

  K_FASTS::getTransfersInter(iptro_tmp, ipt_ndimdx_trans, param_int_tc ,
                                 pair_of_queue);                               

#ifdef TimeShow
  if( rank == 0 )
  {
    time_COM = omp_get_wtime();
    ipt_timecount[4] = ipt_timecount[4] + time_COM -time_init;
    time_init= omp_get_wtime();
  }
#endif

  for (E_Int ircv = 1; ircv < nbcomID +1; ++ircv)
  {
   pt_rcv_queue->emplace_back(param_int_tc[ircv+nbcomIBC + 2], 404);
   CMP::RecvBuffer& recv_buffer = pt_rcv_queue->back_message_buffer();
   recv_buffer.irecv();
  }

#ifdef TimeShow

if( rank == 0 )
  {
    time_COM = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_COM -time_init;  
  }

#endif

  for (E_Int ip2p = 1; ip2p < param_int_tc[0]+1; ++ip2p)
  {
    E_Int ech  = param_int_tc[ip2p+shift_graph];
    dest       = param_int_tc[ech];
    if (dest != rank)  // Inter Process
    {             
      TypeTransfert = 0;
      K_FASTS::setInterpTransfersInter(iptro_tmp,     ipt_ndimdx_trans,  param_int_tc,  param_real_tc ,
                                           param_bci, param_bcf, TypeTransfert , it_target, nidom, ip2p, pair_of_queue, ipt_timecount);   
    }
    // else
    // {
    // TypeTransfert = 0;
    // K_FASTS::setInterpTransfersIntra(iptro_tmp,     ipt_ndimdx_trans,  param_int_tc,  param_real_tc ,
    //                                    param_bci, param_bcf, TypeTransfert , it_target, nidom, ip2p, ipt_timecount); 
    // }
  }
}
// Endif MPI
#endif 

  for (E_Int ip2p = 1; ip2p < param_int_tc[0]+1; ++ip2p)
  {
    E_Int ech  = param_int_tc[ip2p+shift_graph];
    dest       = param_int_tc[ech];
    if (dest == rank)  // Intra Process
    {        
        TypeTransfert = 0;
        K_FASTS::setInterpTransfersIntra(iptro_tmp,     ipt_ndimdx_trans,  param_int_tc,  param_real_tc ,
                                             param_bci, param_bcf, TypeTransfert , it_target, nidom, ip2p, ipt_timecount); 
    }    
  }

// #ifdef TimeShow
//   if(rank == 0)
//   {
//     time_init = omp_get_wtime();
//   }
// #endif

  #ifdef _MPI

  if (init)
  {
      K_FASTS::getTransfersInter(iptro_tmp, ipt_ndimdx_trans, param_int_tc ,
                                     pair_of_queue);

#ifdef TimeShow
  if( rank == 0 )
  {
    time_COM = omp_get_wtime();
    ipt_timecount[4] = ipt_timecount[4] + time_COM -time_init;

    std::cout << "Time in getTransfersInter "     << ipt_timecount[4] << std::endl;
    std::cout << "Time InterpTransfert (Intra)  " << ipt_timecount[1] << std::endl;
    std::cout << "Time in MPI send_buffer, irecv "<< ipt_timecount[0] << std::endl;
    std::cout << "Time InterpTransfert (Inter)  " << ipt_timecount[2] << std::endl;
    std::cout << std::endl << std::endl;
    time_init = omp_get_wtime();
  }    
#endif

      K_FASTS::del_TransferInter(pair_of_queue);
  }
  // MPI Second part (InterCOM ID)
  #endif      

}
}

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
void K_FASTS::setInterpTransfersIntra(
    E_Float**& ipt_ro, E_Int*& ipt_ndimdx, E_Int*& ipt_param_int,
    E_Float*& ipt_param_real, E_Int*& ipt_parambci, E_Float*& ipt_parambcf,
    E_Int& TypeTransfert, E_Int& it_target, E_Int& nidom, E_Int& NoTransfert,
    E_Float*& ipt_timecount)

{
  E_Int rank = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == 0 and timeShowFastS) {
    time_in = omp_get_wtime();
  }

  // Unpack Data
  E_Int varType = ipt_parambci[0];
  E_Int bcType  = ipt_parambci[1];

  E_Float gamma = ipt_parambcf[0];
  E_Float cv = ipt_parambcf[1];
  E_Float muS = ipt_parambcf[2];
  E_Float Cs = ipt_parambcf[3];
  E_Float Ts = ipt_parambcf[4];

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

  E_Int nbcomIBC = ipt_param_int[1];
  E_Int nbcomID  = ipt_param_int[2 + nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 2;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech           = ipt_param_int[NoTransfert + shift_graph];
  E_Int nrac          = ipt_param_int[ech + 1];  // nb total de raccord
  E_Int nrac_inst     = ipt_param_int[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel     = ipt_param_int[ech + 3];  // nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady   = nrac - nrac_inst;        // nb total de raccord stationnaire

//gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0; 
  E_Int pass_inst_fin=1;

  if (nrac_inst > 0) pass_inst_fin=2;


  // on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx = 0;
  for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++) 
  {
    E_Int irac_deb = 0;
    E_Int irac_fin = nrac_steady;
    if(pass_inst == 1){ irac_deb = ipt_param_int[ ech + 4 + it_target             ]; 
                        irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];
                      }

    for (E_Int irac = irac_deb; irac < irac_fin; irac++) 
    {
      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
      if (ipt_param_int[shift_rac + nrac * 10 + 1] > nbRcvPts_mx)
        nbRcvPts_mx = ipt_param_int[shift_rac + nrac * 10 + 1];
    }
  }

  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (bcType <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 13 * threadmax_sdm);
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
          irac_deb = ipt_param_int[ech + 4 + it_target];
          irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
        }

        // printf("iracdeb=  %d, iracfin= %d \n", irac_deb, irac_fin  );
        for (E_Int irac = irac_deb; irac < irac_fin; irac++) {
          // E_Int shift_rac =  ech + 4 + irac;
          E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
          //printf("ipass_typ = %d, pass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,pass_inst,irac , ithread ); 
          // ipass_typ,ipass_inst,irac , ithread );
          E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
          E_Int ibc = 1;
          if (ibcType < 0) ibc = 0;
          if (1 - ibc != ipass_typ) continue;

          E_Int NoD       = ipt_param_int[shift_rac + nrac * 5     ];
          E_Int loc       = ipt_param_int[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
          E_Int NoR       = ipt_param_int[shift_rac + nrac * 11 + 1];
          E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
          E_Int rotation  = ipt_param_int[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal

          E_Int meshtype = 1;  // ONLY FOR STRUCTURE ipt_ndimdxD[NoD +
                               // nidomD*6];
          E_Int cnNfldD = 0;
          E_Int* ptrcnd = NULL;

          for (E_Int eq = 0; eq < nvars_loc; eq++) {
            vectOfRcvFields[eq] = ipt_ro[NoR] + eq * ipt_ndimdx[NoR];
            vectOfDnrFields[eq] = ipt_ro[NoD] + eq * ipt_ndimdx[NoD];
          }
          imd = ipt_ndimdx[NoD + nidom];
          jmd = ipt_ndimdx[NoD + nidom * 2];
          // }

          imdjmd = imd * jmd;

          ////
          //  Interpolation parallele
          ////
          ////

          E_Int nbRcvPts = ipt_param_int[shift_rac + nrac * 10 + 1];
          // E_Int nbDonPts = ipt_param_int[ shift_rac                ];

          E_Int pos;
          pos = ipt_param_int[shift_rac + nrac * 7];      E_Int* ntype    = ipt_param_int  + pos;
          pos = pos + 1 + ntype[0];                       E_Int* types    = ipt_param_int  + pos;
          pos = ipt_param_int[shift_rac + nrac * 6];      E_Int* donorPts = ipt_param_int  + pos;
          pos = ipt_param_int[shift_rac + nrac * 12 + 1]; E_Int* rcvPts   = ipt_param_int  + pos;  // donor et receveur inverser car storage donor
          pos = ipt_param_int[shift_rac + nrac * 8];    E_Float* ptrCoefs = ipt_param_real + pos;

          E_Int nbInterpD = ipt_param_int[shift_rac + nrac];
          E_Float* xPC = NULL;
          E_Float* xPI = NULL;
          E_Float* xPW = NULL;
          E_Float* densPtr = NULL;
          if (ibc == 1) {
            xPC = ptrCoefs + nbInterpD;
            xPI = ptrCoefs + nbInterpD + 3 * nbRcvPts;
            xPW = ptrCoefs + nbInterpD + 6 * nbRcvPts;
            densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;
          }

          E_Int ideb = 0;
          E_Int ifin = 0;
          E_Int shiftCoef = 0;
          E_Int shiftDonnor = 0;

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

                      noi       = shiftDonnor;                             // compteur sur le tableau d indices donneur
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
            #           include "commonInterpTransfers_reorder_neq.h"
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
                                                    densPtr, densPtr+nbRcvPts,
                                                    densPtr+nbRcvPts*2, densPtr+nbRcvPts*3,
                                                    densPtr+nbRcvPts*4, densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, densPtr+nbRcvPts*7, densPtr+nbRcvPts*8,  
                                                    ipt_tmp, size,
                                                    gamma, cv, muS, Cs, Ts,
                                                    vectOfDnrFields, vectOfRcvFields);
                        else if (varType == 2 || varType == 21)
                        {
                          K_CONNECTOR::setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                    xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                                    xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                                    xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                                    densPtr, densPtr+nbRcvPts, 
                                                    densPtr+nbRcvPts*2, densPtr+nbRcvPts*3,
                                                    densPtr+nbRcvPts*4, densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, densPtr+nbRcvPts*7, densPtr+nbRcvPts*8,
                                                    ipt_tmp, size,
                                                    gamma, cv, muS, Cs, Ts,
                                                    vectOfDnrFields, vectOfRcvFields);
                        }
                        else if (varType == 3 || varType == 31)
                          K_CONNECTOR::setIBCTransfersCommonVar3(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                                    xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                                    xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                                    xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                                    densPtr, densPtr +nbRcvPts, 
                                                    densPtr+nbRcvPts*2, densPtr+nbRcvPts*3,
                                                    densPtr+nbRcvPts*4, densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, densPtr+nbRcvPts*7, densPtr+nbRcvPts*8,
                                                    ipt_tmp, size,
                                                    gamma, cv, muS, Cs, Ts,
                                                    vectOfDnrFields, vectOfRcvFields);
                      }//ibc          
              //*
              //        } //chunk
              //*/
                      ideb       = ideb + ifin;
                      shiftCoef  = shiftCoef   +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
                      shiftDonnor= shiftDonnor +  ntype[1+ndtyp];           //shift donnor entre 2 types successif
                   }// type 
                }//irac
               }//pass_inst
              #pragma omp barrier 
    }  // ipass
  }    // omp

  if (rank == 0 and timeShowFastS) {
    time_out = omp_get_wtime();
    ipt_timecount[1] = ipt_timecount[1] + time_out - time_in;
    time_in = omp_get_wtime();
  }

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
    E_Float**& ipt_ro, E_Int*& ipt_ndimdx, E_Int*& ipt_param_int,
    E_Float*& ipt_param_real, E_Int*& ipt_parambci, E_Float*& ipt_parambcf,
    E_Int& TypeTransfert, E_Int& it_target, E_Int& nidom, E_Int& NoTransfert,
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue, E_Float*& ipt_timecount)

{
  E_Int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0 and timeShowFastS) {
    time_in = omp_get_wtime();
  }

  // Unpack Data
  E_Int varType = ipt_parambci[0];
  E_Int bcType = ipt_parambci[1];

  E_Float gamma = ipt_parambcf[0];
  E_Float cv = ipt_parambcf[1];
  E_Float muS = ipt_parambcf[2];
  E_Float Cs = ipt_parambcf[3];
  E_Float Ts = ipt_parambcf[4];

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

  E_Int nbcomIBC = ipt_param_int[1];
  E_Int nbcomID  = ipt_param_int[2 + nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 2;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech  = ipt_param_int[NoTransfert + shift_graph];
  E_Int nrac = ipt_param_int[ech + 1];  // nb total de raccord
  E_Int nrac_inst = ipt_param_int[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel = ipt_param_int[ech + 3];  // nb de pas de temps stocker pour
                                             // chaque raccord instationnaire

  
  E_Int nrac_steady = nrac - nrac_inst;  // nb total de raccord stationnaire
  E_Int pass_inst_deb=0; 
  E_Int pass_inst_fin=1;
  E_Int nrac_inst_level = 0;
  if (rank == 0 and timeShowFastS ) {
    time_in = omp_get_wtime();
  }

  if (nrac_inst > 0) {
   pass_inst_fin=2;
   nrac_inst_level = ipt_param_int[ech + 4 + it_target + timelevel] - ipt_param_int[ech + 4 + it_target] + 1; 
  } 
  
  // on dimension tableau travail pour IBC et pour transfert
  // E_Int nrac_inst_level = ipt_param_int[ech + 4 + it_target + timelevel] -
  //                         ipt_param_int[ech + 4 + it_target] + 1;

  std::vector<E_Float*> frp(nrac_steady + nrac_inst_level);

  int dest = ipt_param_int[ipt_param_int[NoTransfert + shift_graph]]; 

  SendQueue* pt_snd_queue = pair_of_queue->second;
  pt_snd_queue->emplace_back(dest, 404);

  // Préparation du buffer d'envoi :
  CMP::SendBuffer& send_buffer = pt_snd_queue->back_message_buffer();
  // A partir d'ici pour allouer les tableaux à remplir
  std::vector<CMP::SendBuffer::PackedData*> pck_data;

  
  // Test IF Has Data to send and count rac:
  // Nb rac.
  bool has_data_to_send = false;
  E_Int count_rac = 0;
  E_Int nbRcvPts_mx = 0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++) 
    {
        E_Int irac_deb = 0;
        E_Int irac_fin = nrac_steady;
        if ( pass_inst == 1 ) 
        {
            irac_deb = ipt_param_int[ech + 4 + it_target];
            irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
        }
        for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
        {
  // for (E_Int irac = 0; irac < nrac; irac++) {
        E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
        E_Int nbRcvPts = ipt_param_int[shift_rac + nrac * 10 + 1];
    
        E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
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
        }
    }   

  if (has_data_to_send) {
    if (rank == 0 and timeShowFastS) {
    time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
    time_in  = omp_get_wtime();
  }
    
    send_buffer << count_rac;

    count_rac = 0;
    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) {
          irac_deb = ipt_param_int[ech + 4 + it_target];
          irac_fin = ipt_param_int[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
      {
    
      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
      
      E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
      E_Int ibc = 1;
      if (ibcType < 0) ibc = 0;

      if (TypeTransfert == 0 && ibc == 1) {
        continue;
      } else if (TypeTransfert == 1 && ibc == 0) {
        continue;
      }
      E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];
      E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];  // flag raccord rans/LES
      E_Int Nozone = ipt_param_int[shift_rac + nrac * 11 + 1];

      send_buffer << Nozone;  

      pck_data.push_back(&send_buffer.push_inplace_array(nvars_loc * nbRcvPts *
                                                         sizeof(E_Float)));    
      E_Int PtlistDonor  = ipt_param_int[shift_rac + nrac * 12 + 1];
      E_Int* ipt_listRcv = ipt_param_int + PtlistDonor;

      send_buffer << CMP::vector_view<E_Int>(ipt_listRcv, nbRcvPts);
      count_rac += 1;
      }
    }
    
    send_buffer.finalize_and_copy();

    count_rac = 0;
    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) {
          irac_deb = ipt_param_int[ech + 4 + it_target];
          irac_fin = ipt_param_int[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
      {
      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
      E_Int ibcType   = ipt_param_int[shift_rac + nrac * 3];
      E_Int ibc = 1;
      if (ibcType < 0) ibc = 0;

      if (TypeTransfert == 0 && ibc == 1) {
        continue;
      } else if (TypeTransfert == 1 && ibc == 0) {
        continue;
      }

      frp[count_rac] = pck_data[count_rac]->data<E_Float>();
      count_rac += 1;
      }
    }    
  }


  if (rank == 0 and timeShowFastS) {
    time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
    time_in  = omp_get_wtime();
  }

  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (bcType <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 13 * threadmax_sdm);
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
    E_Int indD0, indD, i, j, k, ncfLoc, nocf, indCoef, noi, sizecoefs, Nbchunk,
        imd, jmd, imdjmd;

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    // 1ere pass: IBC
    // 2eme pass: transfert
    //
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
        irac_deb = ipt_param_int[ech + 4 + it_target];
                    irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
                }

                for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) 
                {
                   E_Int shift_rac = ech + 4 + timelevel*2 + irac;
                   E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
                   E_Int ibc = 1;
                   if (ibcType < 0) ibc = 0;
                   if (1 - ibc != ipass_typ) continue;
         
                   E_Int NoD       = ipt_param_int[shift_rac + nrac * 5     ];
                   E_Int loc       = ipt_param_int[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
                   E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];
                   E_Int NoR       = ipt_param_int[shift_rac + nrac * 11 + 1];
                   E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
                   E_Int rotation  = ipt_param_int[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal
         
                   E_Int meshtype = 1;  // ONLY FOR STRUCTURE ipt_ndimdxD[NoD + nidom*6];
                   E_Int cnNfldD  = 0;
                   E_Int* ptrcnd  = NULL;
         
                   // printf("navr_loc %d %d %d \n", nvars_loc, nvars, Rans);
                   for (E_Int eq = 0; eq < nvars_loc; eq++) {
                     vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                     vectOfDnrFields[eq] = ipt_ro[NoD] + eq * ipt_ndimdx[NoD];
                   }
                   imd = ipt_ndimdx[NoD + nidom];
                   jmd = ipt_ndimdx[NoD + nidom * 2];                  

                   imdjmd = imd * jmd;
                   ////
                   //  Interpolation parallele
                   ////
                   ////
                   E_Int pos;
                   pos               = ipt_param_int[shift_rac + nrac * 7];
                   E_Int* ntype      = ipt_param_int + pos;
                   pos               = pos + 1 + ntype[0];
                   E_Int* types      = ipt_param_int + pos;
                   pos               = ipt_param_int[shift_rac + nrac * 6];
                   E_Int* donorPts   = ipt_param_int + pos;
                   pos               = ipt_param_int[shift_rac + nrac * 8];
                   E_Float* ptrCoefs = ipt_param_real + pos;
                   // pos               = ipt_param_int[shift_rac + nrac * 12 + 1];
                   // E_Int* rcvPts     = ipt_param_int +  pos;
         
                   E_Int nbInterpD = ipt_param_int[shift_rac + nrac];
                   E_Float* xPC = NULL;
                   E_Float* xPI = NULL;
                   E_Float* xPW = NULL;
                   E_Float* densPtr = NULL;
                   if (ibc == 1) {
                     xPC = ptrCoefs + nbInterpD;
                     xPI = ptrCoefs + nbInterpD + 3 * nbRcvPts;
                     xPW = ptrCoefs + nbInterpD + 6 * nbRcvPts;
                     densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;
                   }
         
                   E_Int ideb = 0;
                   E_Int ifin = 0;
                   E_Int shiftCoef = 0;
                   E_Int shiftDonnor = 0;
         
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
         
         
                  noi     = shiftDonnor;  // compteur sur le tableau d indices donneur
                  indCoef = ( pt_deb - ideb ) * sizecoefs + shiftCoef;
                  E_Int NoR = ipt_param_int[shift_rac + nrac * 11 + 1];
                  //if (ipt_param_int[ech]==0) printf("No rac= %d , NoR= %d, NoD= %d, Ntype= %d, ptdeb= %d, ptfin= %d, NptD= %d, neq= %d, skip= %d, rank= %d, dest= %d,  thread= %d\n",
                  //irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin , 
                  //ipt_param_int[ shift_rac + nrac*10+1  ], ipt_param_int[ shift_rac + nrac*13+1  ], ipt_param_int[ shift_rac + nrac*15+1  ], 
                  //rank, ipt_param_int[ ech  ], ithread );
                  if ( nvars_loc == 5 ) {
#include "commonInterpTransfersD_reorder_5eq.h"
                        } else if ( nvars_loc == 6 ) {
#include "commonInterpTransfersD_reorder_6eq.h"
                        } else {
#include "commonInterpTransfersD_reorder_neq.h"
                        }
                  // Prise en compte de la periodicite par rotation
                  if ( rotation == 1 ) {
                      E_Float* angle = ptrCoefs + nbInterpD;
#include "includeTransfersD_rotation.h"
                        }
                  // ibc
                  if (ibc == 1) {
                      // tableau temporaire pour utiliser la routine commune K_CONNECTOR::setIBCTransfersCommon
                      for ( E_Int noind = pt_deb; noind < pt_fin; noind++ ) rcvPts[noind] = noind;
                      if ( varType == 1 || varType == 11 )
                          K_CONNECTOR::setIBCTransfersCommonVar1( bcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, xPC,
                                                     xPC + nbRcvPts, xPC + nbRcvPts * 2, xPW, xPW + nbRcvPts,
                                                     xPW + nbRcvPts * 2, xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
                                                     densPtr, densPtr+nbRcvPts, 
                                                     densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, 
                                                     densPtr+nbRcvPts*4, densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, densPtr+nbRcvPts*7, densPtr+nbRcvPts*8,
                                                     ipt_tmp, size, gamma, cv, muS, Cs,
                                                     Ts, vectOfDnrFields, vectOfRcvFields );
                      else if ( varType == 2 || varType == 21 )
                          K_CONNECTOR::setIBCTransfersCommonVar2( bcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, xPC,
                                                     xPC + nbRcvPts, xPC + nbRcvPts * 2, xPW, xPW + nbRcvPts,
                                                     xPW + nbRcvPts * 2, xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
                                                     densPtr, densPtr+nbRcvPts, 
                                                     densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, 
                                                     densPtr+nbRcvPts*4, densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, densPtr+nbRcvPts*7, densPtr+nbRcvPts*8,  
                                                     ipt_tmp, size, gamma, cv, muS, Cs,
                                                     Ts, vectOfDnrFields, vectOfRcvFields );
                      else if ( varType == 3 || varType == 31 )
                          K_CONNECTOR::setIBCTransfersCommonVar3( bcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, xPC,
                                                     xPC + nbRcvPts, xPC + nbRcvPts * 2, xPW, xPW + nbRcvPts,
                                                     xPW + nbRcvPts * 2, xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
                                                     densPtr, densPtr+nbRcvPts, 
                                                     densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, 
                                                     densPtr+nbRcvPts*4, densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, densPtr+nbRcvPts*7, densPtr+nbRcvPts*8,  
                                                     ipt_tmp, size, gamma, cv, muS, Cs,
                                                     Ts, vectOfDnrFields, vectOfRcvFields );
                  }  // ibc
                  //        } //chunk
                  ideb        = ideb + ifin;
                  shiftCoef   = shiftCoef + ntype[1 + ndtyp] * sizecoefs;  // shift coef   entre 2 types successif
                  shiftDonnor = shiftDonnor + ntype[1 + ndtyp];            // shift donnor entre 2 types successif
                  }                                                        // type
                    count_rac += 1;
                }  // irac
            }      // pass_inst
#pragma omp barrier
    }  // ipass
  }    // omp

  if (rank == 0 and timeShowFastS) {
    time_out = omp_get_wtime();
    ipt_timecount[2] = ipt_timecount[2] + time_out - time_in;
    time_in = omp_get_wtime();
  }

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
  
  if (rank == 0 and timeShowFastS) {
    time_out = omp_get_wtime();
    ipt_timecount[0] = ipt_timecount[0] + time_out - time_in;
  }



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
void K_FASTS::getTransfersInter(
    E_Float**& ipt_roD, E_Int*& ipt_ndimdxD, E_Int*& ipt_param_int,
    std::pair<RecvQueue*, SendQueue*>*& pair_of_queue) {
  
  E_Int rank;

  #ifdef _MPI
  int init = 0;
  if (init) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  #endif

  // Attente finalisation de la réception :
  assert(pair_of_queue != NULL);

  vector<CMP::vector_view<E_Float> > recv_frp(2048);
  vector<E_Int> recv_nozone(2048);
  vector<E_Int> recv_nvarloc(2048);
  vector<CMP::vector_view<E_Int> > recv_listRc(2048);

  RecvQueue* pt_rcv_queue = pair_of_queue->first;
  SendQueue* pt_snd_queue = pair_of_queue->second;

  while (not pt_rcv_queue->empty()) {

    RecvQueue::iterator it = pt_rcv_queue->get_first_complete_message();
    if (it != pt_rcv_queue->end()) {  // ok, une réception de prête*/

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
          E_Int decal = ipt_ndimdxD[recv_nozone[irac]];

          if (recv_nvarloc[irac] == 5) {
#pragma omp for
            for (int irecv = 0; irecv < sz; ++irecv) {
              //      for (int eq=0; eq < recv_nvarloc[irac] ; ++eq){
              ilistrecv =
                  recv_listRc[irac]
                             [irecv];  //+ eq*ipt_ndimdxD[recv_nozone[irac]];
              ipt_roD[recv_nozone[irac]][ilistrecv] = recv_frp[irac][irecv];
              ipt_roD[recv_nozone[irac]][ilistrecv + decal] =
                  recv_frp[irac][irecv + 1 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 2 * decal] =
                  recv_frp[irac][irecv + 2 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 3 * decal] =
                  recv_frp[irac][irecv + 3 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 4 * decal] =
                  recv_frp[irac][irecv + 4 * sz];
              //}
            }  // end for (int irecv
          } else if (recv_nvarloc[irac] == 6) {
#pragma omp for
            for (int irecv = 0; irecv < sz; ++irecv) {
              ilistrecv =
                  recv_listRc[irac]
                             [irecv];  //+ eq*ipt_ndimdxD[recv_nozone[irac]];
              ipt_roD[recv_nozone[irac]][ilistrecv] = recv_frp[irac][irecv];
              ipt_roD[recv_nozone[irac]][ilistrecv + decal] =
                  recv_frp[irac][irecv + 1 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 2 * decal] =
                  recv_frp[irac][irecv + 2 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 3 * decal] =
                  recv_frp[irac][irecv + 3 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 4 * decal] =
                  recv_frp[irac][irecv + 4 * sz];
              ipt_roD[recv_nozone[irac]][ilistrecv + 5 * decal] =
                  recv_frp[irac][irecv + 5 * sz];
            }  // end for (int irecv
          }

        }  // end omp parallel
      }  // end for ( int irac = ....

#if defined(PROFILE_TRANSFERT)
      end = omp_get_wtime();
      recv_buffer_time += end - beg;
#endif
      // delete [] ipt_ndimdxD; delete [] ipt_roD;
      pt_rcv_queue->pop(it);
    }  // end if (it != pt_msg_manager->end()infos
  }    // End  while (not pt_msg_manager->empty())

  pt_snd_queue->waitAll();
  pt_snd_queue->clear();
#if defined(PROFILE_TRANSFERT)
  std::cout << "Temps actuel passe :  " << send_buffer_time << "( " << isend_time
            << " ) vs " << recv_buffer_time << std::endl;
#endif
}
#endif