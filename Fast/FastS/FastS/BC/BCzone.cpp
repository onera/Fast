/*
    Copyright 2013-2025 Onera.

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
#include "FastS/fastS.h"
#include "FastS/param_solver.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

// default state: staref de la base ou zone  si pas de staref dans la bc
#define DEFAULT_STATE(string)                                       \
  FldArrayF data(6);                                                \
  E_Float *ipt_data_loc = data.begin();                             \
  if (nbdata == 0)                                                  \
  {                                                                 \
    ipt_data_loc[0] = param_real[ROINF];                            \
    ipt_data_loc[1] = param_real[VXINF] * param_real[ROINF];        \
    ipt_data_loc[2] = param_real[VYINF] * param_real[ROINF];        \
    ipt_data_loc[3] = param_real[VZINF] * param_real[ROINF];        \
    ipt_data_loc[4] = param_real[PINF] / (param_real[GAMMA] - 1.) + \
                      0.5 * (param_real[VXINF] * ipt_data_loc[1] +  \
                             param_real[VYINF] * ipt_data_loc[2] +  \
                             param_real[VZINF] * ipt_data_loc[3]);  \
    ipt_data_loc[5] = param_real[RONUTILDEINF];                     \
  }                                                                 \
  else if (nbdata == 21)                                            \
  {                                                                 \
    ipt_data_loc[0] = ipt_data[3];                                  \
    ipt_data_loc[1] = ipt_data[4];                                  \
    ipt_data_loc[2] = ipt_data[5];                                  \
    ipt_data_loc[3] = ipt_data[6];                                  \
    ipt_data_loc[4] = ipt_data[7];                                  \
    ipt_data_loc[5] = ipt_data[16];                                 \
  }                                                                 \
  else                                                              \
  {                                                                 \
    ipt_data_loc = ipt_data;                                        \
  }

using namespace std;
using namespace K_FLD;

//
//=============================================================================
E_Int K_FASTS::BCzone(
    E_Int &nd, E_Int &lrhs, E_Int &nstep, E_Int &lcorner,
    E_Int *param_int, E_Float *param_real, E_Int &npass, E_Int *ipt_ind_dm, E_Int *ipt_ind_dm_thread, E_Int *ipt_ind_CL,
    E_Int *ipt_ind_CL119, E_Int *ipt_ind_CLgmres, E_Int *ishift_lu, E_Float *iptrop, E_Float *ipti, E_Float *iptj, E_Float *iptk, E_Float *iptx,
    E_Float *ipty, E_Float *iptz, E_Float *iptventi, E_Float *iptventj, E_Float *iptventk, E_Float *iptrop_gmres, E_Float *iptmut)
{

  E_Int err = 1;
  E_Int neq_mtr;
  E_Int itest[6], lskip[6], ind_rhs[6], ind_avg[6], ind_avg_thread[6], ind_mjr[6], ind_mjr_thread[6];
  E_Int lrhs_loc = 0;
  if (lrhs == 1)
    lrhs_loc = 1;

  itest[0] = 1;
  itest[1] = param_int[IJKV];
  itest[2] = 1;
  itest[3] = param_int[IJKV + 1];
  itest[4] = 1;
  itest[5] = param_int[IJKV + 2];

  vector<PyArrayObject *> hook;

  for (E_Int i = 0; i < 6; i++)
    lskip[i] = 1;

  E_Int *ipt_nijk = param_int + NIJK;
  E_Int *ipt_ijkv = param_int + IJKV;

  E_Float *iptijk;
  E_Float *ipventijk;

  E_Int pt_bcs = param_int[PT_BC];
  E_Int nb_bc = param_int[pt_bcs];

  E_Int ificmax[6];
  if (param_int[LU_MATCH] != 0)
    param_int[LU_MATCH] = 1;

  for (E_Int ipara = 0; ipara < 4; ipara++)
  {
    ificmax[ipara] = param_int[NIJK + 3] * param_int[LU_MATCH];
  }
  for (E_Int ipara = 4; ipara < 6; ipara++)
  {
    ificmax[ipara] = param_int[NIJK + 4] * param_int[LU_MATCH];
  }

  for (E_Int ndf = 0; ndf < nb_bc; ndf++)
  {
    E_Int pt_bc = param_int[pt_bcs + 1 + ndf];
    E_Int idir = param_int[pt_bc + BC_IDIR];

    E_Int ipara = idir - 1;

    // on test si le sous domaine touche le bord du domaine et si la CL est bonne candidate a implicitation
    if (lrhs_loc == 1)
    {
      if ((ipt_ind_dm_thread[ipara] == itest[ipara]) &&
          ((param_int[pt_bc + BC_TYPE] >= 3 && param_int[pt_bc + BC_TYPE] <= 7) || param_int[pt_bc + BC_TYPE] == 12 || param_int[pt_bc + BC_TYPE] == 18 || param_int[pt_bc + BC_TYPE] == BCWALLMODEL
           //||   param_int[pt_bc + BC_TYPE] == 31
           ))
      {
        lskip[ipara] = 0;
        lcorner = 0;
        if ((param_int[pt_bc + BC_TYPE] >= 3 && param_int[pt_bc + BC_TYPE] <= 7) || param_int[pt_bc + BC_TYPE] == 12)
          ificmax[ipara] = 1;
      }
    }
    else
    {
      lskip[ipara] = 1;
    }

  } // fin loop CL warmup

  if (lrhs_loc == 1)
  {
    ishift_lu[0] = ipt_ind_dm_thread[0] - (1 - lskip[0]) * ificmax[0];
    ishift_lu[1] = ipt_ind_dm_thread[1] + (1 - lskip[1]) * ificmax[1];
    ishift_lu[2] = ipt_ind_dm_thread[2] - (1 - lskip[2]) * ificmax[2];
    ishift_lu[3] = ipt_ind_dm_thread[3] + (1 - lskip[3]) * ificmax[3];
    ishift_lu[4] = ipt_ind_dm_thread[4] - (1 - lskip[4]) * ificmax[4];
    ishift_lu[5] = ipt_ind_dm_thread[5] + (1 - lskip[5]) * ificmax[5];

    if (param_int[LU_MATCH] == 1)
    {
      E_Int nfic_ij = param_int[NIJK + 3];

      if (ishift_lu[0] > ipt_ind_dm[0])
      {
        ishift_lu[0] -= nfic_ij;
      }
      if (ishift_lu[1] < ipt_ind_dm[1])
      {
        ishift_lu[1] += nfic_ij;
      }
      if (ishift_lu[2] > ipt_ind_dm[2])
      {
        ishift_lu[2] -= nfic_ij;
      }
      if (ishift_lu[3] < ipt_ind_dm[3])
      {
        ishift_lu[3] += nfic_ij;
      }
      // if (ishift_lu[0]> 1                  ) {ishift_lu[0] -= nfic_ij;}
      // if (ishift_lu[1]< param_int[ IJKV   ]) {ishift_lu[1] += nfic_ij;}
      // if (ishift_lu[2]> 1                  ) {ishift_lu[2] -= nfic_ij;}
      // if (ishift_lu[3]< param_int[ IJKV+1 ]) {ishift_lu[3] += nfic_ij;}
      if (param_int[ITYPZONE] != 3)
      {
        E_Int nfic_k = param_int[NIJK + 4];

        if (ishift_lu[4] > ipt_ind_dm[4])
        {
          ishift_lu[4] -= nfic_k;
        }
        if (ishift_lu[5] < ipt_ind_dm[5])
        {
          ishift_lu[5] += nfic_k;
        }
        // if (ishift_lu[4] > 1                    ){ ishift_lu[4] -= nfic_k;}
        // if (ishift_lu[5] < param_int[ IJKV + 2 ]){ ishift_lu[5] += nfic_k;}
      }
    }

    E_Int idirmax = 6;
    if (param_int[ITYPZONE] == 3)
      idirmax = 4;

    // extrapolation du RHS par defaut sur toutes les faces ou on etend la resolution LU.  (Sinon probleme si face
    // domaine= paroi + raccord)
    for (E_Int idir = 0; idir < idirmax; idir++)
    {
      for (E_Int l = 0; l < 6; l++)
        ind_rhs[l] = ipt_ind_dm_thread[l];
      E_Int lskip_loc;
      if (idir == 0)
      {
        ind_rhs[0] = 1;
        ind_rhs[1] = 1;
      }
      if (idir == 2)
      {
        ind_rhs[2] = 1;
        ind_rhs[3] = 1;
      }
      if (idir == 4)
      {
        ind_rhs[4] = 1;
        ind_rhs[5] = 1;
      }
      if (idir == 1)
      {
        ind_rhs[0] = ind_rhs[1];
      }
      if (idir == 3)
      {
        ind_rhs[2] = ind_rhs[3];
      }
      if (idir == 5)
      {
        ind_rhs[4] = ind_rhs[5];
      }

      if (lskip[idir] == 0)
      {
        E_Int idir_loc = idir + 1;
        E_Int typebc_loc = 0;
        indice_cl_sdm_(idir_loc, npass, lskip_loc, typebc_loc, ipt_nijk[3], ipt_nijk[4], // IN
                       ipt_ijkv, ind_rhs, ipt_ind_dm, ipt_ind_dm_thread,                 // IN
                       ipt_ind_CL, ipt_ind_CL119);

        if (lskip_loc == 0)
        {
          E_Int nbdata = 0;
          E_Float mobile_coef = 1.;
          E_Float *ipt_data;
          DEFAULT_STATE("BCExtrapolate")
          // MUSCL
          E_Float c4, c5, c6;
          c4 = 5. / 6.;
          c5 = 2. / 6.;
          c6 = -1. / 6.;

          bvbs_extrapolate_(idir_loc, lrhs_loc, nstep,
                            neq_mtr, mobile_coef,
                            param_int, ipt_ind_CL119,
                            param_real, c4, c5, c6,
                            iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
      }
    }
  }

  for (E_Int ndf = 0; ndf < nb_bc; ndf++)
  {
    E_Int pt_bc = param_int[pt_bcs + 1 + ndf];

    E_Int idir = param_int[pt_bc + BC_IDIR];
    E_Int nbdata = param_int[pt_bc + BC_NBDATA];
    E_Int bc_type = param_int[pt_bc + BC_TYPE];

    // printf("ptbc= %d , idir= %d , ndom= %d , bctype=  %d, nstep= %d , nbdata= %d , lrhs= %d \n", pt_bc, idir,nd,  bc_type, nstep,nbdata, lrhs_loc );

    E_Int *iptsize_data = param_int + pt_bc + BC_NBDATA + 1;

    E_Int ipara = idir - 1;

    if (lskip[ipara] == 0 || lrhs_loc == 0) // test pour savoir si face ipara est skipper en LU pour le rhs
    {
      if (idir <= 2)
      {
        iptijk = ipti;
        neq_mtr = param_int[NEQ_IJ];
        ipventijk = iptventi;
      }
      else if (idir <= 4)
      {
        iptijk = iptj;
        neq_mtr = param_int[NEQ_IJ];
        ipventijk = iptventj;
      }
      else
      {
        iptijk = iptk;
        neq_mtr = param_int[NEQ_K];
        ipventijk = iptventk;
      }

      // calcul intersection sous domaine et fenetre CL
      E_Int lskip_loc;
      indice_cl_sdm_(idir, npass, lskip_loc, param_int[pt_bc + BC_TYPE], ipt_nijk[3], ipt_nijk[4], // IN
                     ipt_ijkv, param_int + pt_bc + BC_FEN, ipt_ind_dm, ipt_ind_dm_thread,          // IN
                     ipt_ind_CL, ipt_ind_CL119);                                                   // OUT

      E_Float *ipt_data;

      if (lrhs == 2)
      {
        E_Int range = 2;
        E_Int without_ghost = 1;
        E_Float signe = 1.;
        // E_Int ipt_ind_CLgmres[6];

        for (E_Int i = 0; i < 6; i++)
          ipt_ind_CLgmres[i] = ipt_ind_CL[i];

        if (idir % 2 == 1)
        {
          ipt_ind_CLgmres[idir - 1] += range * without_ghost;
          ipt_ind_CLgmres[idir] += range;
        }
        else
        {
          // printf("idir2   %d \n", idir);
          ipt_ind_CLgmres[idir - 2] -= range;
          ipt_ind_CLgmres[idir - 1] -= range * without_ghost;
        }
        // printf("idir  %d %d %d %d  %d %d %d\n", idir,  ipt_ind_CLgmres[0], ipt_ind_CLgmres[1], ipt_ind_CLgmres[2],  ipt_ind_CLgmres[3], ipt_ind_CLgmres[4], ipt_ind_CLgmres[5]);

        pre_bc_(param_int, signe, ipt_ind_CLgmres, iptrop, iptrop_gmres);
      }

      if (nbdata != 0)
        ipt_data = param_real + param_int[pt_bcs + 1 + ndf + nb_bc];

      if (lskip_loc == 0)
      {

        if (lrhs_loc == 1 && (bc_type == 0 || bc_type == 1 || bc_type == 2 || bc_type == 10 || bc_type == 11 ||
                              bc_type == 14 || bc_type == 16 || bc_type == 17 || bc_type == 20 || bc_type == 21))
          bc_type = 0; // si rhs, on extrapole sauf si wall

        // MUSCL
        E_Float c4, c5, c6;
        c4 = 5. / 6.;
        c5 = 2. / 6.;
        c6 = -1. / 6.;

        if (bc_type == 0)
        {
          DEFAULT_STATE("BCExtrapolate")
          E_Float mobile_coef = 1.;
          bvbs_extrapolate_(idir, lrhs_loc, nstep,
                            neq_mtr, mobile_coef,
                            param_int, ipt_ind_CL119,
                            param_real, c4, c5, c6,
                            iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 1)
        {
          DEFAULT_STATE("BCFarfield")
          E_Float mobile_coef = 1.;

          bvbs_farfield_(idir, lrhs_loc, nstep,
                         neq_mtr, mobile_coef,
                         param_int, ipt_ind_CL,
                         param_real, c4, c5, c6,
                         iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 10)
        {
          DEFAULT_STATE("BCOutflow")
          E_Float mobile_coef = 1.;

          bvbs_outflow_(idir, lrhs_loc, nstep,
                        neq_mtr, mobile_coef,
                        param_int, ipt_ind_CL,
                        param_real, c4, c5, c6,
                        iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 13)
        {
          DEFAULT_STATE("BCInflow")
          E_Float mobile_coef = 1.;

          if (nbdata == 0 || (nbdata != 0 && iptsize_data[0] <= 6))
          {
            bvbs_inflow_(idir, lrhs_loc, nstep,
                         neq_mtr, mobile_coef,
                         param_int, ipt_ind_CL,
                         param_real, c4, c5, c6,
                         iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
          }
          else
          {
#include "BC/INCREMENT_BC.h"
            E_Float mobile_coef = 1.;

            if (nbdata != 6 and param_int[IFLOW] == 3)
            {
              printf("ERROR BCInflow: inappropriate number of BC data \n");
              exit(err);
            }
            if (nbdata != 5 and param_int[IFLOW] <= 2)
            {
              printf("ERROR BCInflow: inappropriate number of BC data \n");
              exit(err);
            }

            E_Float *ipt_data1 = ipt_data;
            E_Float *ipt_data2 = ipt_data1 + iptsize_data[0];
            E_Float *ipt_data3 = ipt_data2 + iptsize_data[0];
            E_Float *ipt_data4 = ipt_data3 + iptsize_data[0];
            E_Float *ipt_data5 = ipt_data4 + iptsize_data[0];
            E_Float *ipt_data6 = ipt_data5 + iptsize_data[0];

            bvbs_inflow_fich_(idir, lrhs_loc, nstep,
                              neq_mtr, mobile_coef,
                              param_int, ipt_ind_CL,
                              param_real, c4, c5, c6,
                              iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                              ipt_data1, ipt_data2, ipt_data3, ipt_data4, ipt_data5, ipt_data6,
                              iptsize_data[0], inc_bc);
          }
        }
        else if (bc_type == 14)
        { // extrap nutilde interface RANSLES
          DEFAULT_STATE("BCExtrapolate")
          E_Float mobile_coef = 1.;
          bvbs_extrapolate_ransles_(idir, lrhs_loc, nstep,
                                    neq_mtr, mobile_coef,
                                    param_int, ipt_ind_CL119,
                                    param_real, c4, c5, c6,
                                    iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 19)
        { // inflow Lund
          DEFAULT_STATE("BCInflowLund")

#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1.;

          E_Float *ipt_in1 = ipt_data;
          E_Float *ipt_in2 = ipt_in1 + iptsize_data[0];
          E_Float *ipt_in3 = ipt_in2 + iptsize_data[0];
          E_Float *ipt_in4 = ipt_in3 + iptsize_data[0];
          E_Float *ipt_in5 = ipt_in4 + iptsize_data[0];

          E_Float *ipt_pl1 = ipt_data + iptsize_data[0] * param_int[NEQ];
          E_Float *ipt_pl2 = ipt_pl1 + iptsize_data[0];
          E_Float *ipt_pl3 = ipt_pl2 + iptsize_data[0];
          E_Float *ipt_pl4 = ipt_pl3 + iptsize_data[0];
          E_Float *ipt_pl5 = ipt_pl4 + iptsize_data[0];

          E_Float *ipt_pl6 = ipt_pl5;
          E_Float *ipt_in6 = ipt_in5;
          if (param_int[NEQ] == 6)
          {
            ipt_pl6 = ipt_pl5 + iptsize_data[0];
            ipt_in6 = ipt_in5 + iptsize_data[0];
          }

          E_Float *ipt_paramlund = ipt_data + 2 * iptsize_data[0] * param_int[NEQ];

          bvbs_inflow_lund_(idir, lrhs_loc, nstep,
                            neq_mtr, mobile_coef,
                            param_int, ipt_ind_CL,
                            param_real, c4, c5, c6,
                            iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                            ipt_in1, ipt_in2, ipt_in3, ipt_in4, ipt_in5, ipt_in6,
                            ipt_pl1, ipt_pl2, ipt_pl3, ipt_pl4, ipt_pl5, ipt_pl6, ipt_paramlund, iptsize_data[0], inc_bc);

          // Added romain Paris           /////////////////////////////////////////////////////
        }
        else if (bc_type == 20)
        { // inj MFR Elsa-like injection massflow rate BC
          DEFAULT_STATE("BCInjMFR")

#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1.;

          if (nbdata != 6 and param_int[IFLOW] == 3)
          {
            printf("ERROR BCInjMFR: inappropriate number of BC data \n");
            exit(err);
          }
          if (nbdata != 5 and param_int[IFLOW] <= 2)
          {
            printf("ERROR BCInjMFR: inappropriate number of BC data \n");
            exit(err);
          }

          E_Float *ipt_d0x = ipt_data;                       // d0x
          E_Float *ipt_d0y = ipt_data + iptsize_data[0];     // d0y
          E_Float *ipt_d0z = ipt_data + 2 * iptsize_data[0]; // d0z
          E_Float *ipt_qp = ipt_data + 3 * iptsize_data[0];  // qp
          E_Float *ipt_ha = ipt_data + 4 * iptsize_data[0];  // ha
          E_Float *ipt_nue = ipt_data + 5 * iptsize_data[0]; // nue

          bvbs_injmfr_(idir, lrhs_loc, nstep,
                       neq_mtr, mobile_coef,
                       param_int, ipt_ind_CL,
                       param_real, c4, c5, c6,
                       iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                       ipt_d0x, ipt_d0y, ipt_d0z, ipt_qp, ipt_ha, ipt_nue, iptsize_data[0], inc_bc);
        }
        else if (bc_type == 21)
        { // out MFR
          DEFAULT_STATE("BCOutMFR")

#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1.;

          if (nbdata < 1)
          {
            printf("ERROR BCOutMFR: inappropriate number of BC data \n");
            exit(err);
          }

          E_Float *ipt_qp = ipt_data;

          bvbs_outmfr_(idir, lrhs_loc, nstep,
                       neq_mtr, mobile_coef,
                       param_int, ipt_ind_CL,
                       param_real, c4, c5, c6,
                       iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                       ipt_qp, iptsize_data[0], inc_bc);
        }
        else if (bc_type == 2)
        {
          DEFAULT_STATE("BCSupersonic")
          E_Float mobile_coef = 1.;

          if (nbdata == 0 || (nbdata != 0 && iptsize_data[0] <= 6))
          {
            bvbs_inflow_supersonic_(idir, lrhs_loc, nstep,
                                    neq_mtr, mobile_coef,
                                    param_int, ipt_ind_CL,
                                    param_real, c4, c5, c6,
                                    iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
          }
          else
          {
#include "BC/INCREMENT_BC.h"
            if (nbdata != 6 and param_int[IFLOW] == 3)
            {
              printf("ERROR BCInflowSupersonic: inappropriate number of BC data \n");
              exit(err);
            }
            if (nbdata != 5 and param_int[IFLOW] <= 2)
            {
              printf("ERROR BCInflowSupersonic: inappropriate number of BC data \n");
              exit(err);
            }

            E_Float *ipt_data1 = ipt_data;                    // ro
            E_Float *ipt_data2 = ipt_data1 + iptsize_data[0]; // rou
            E_Float *ipt_data3 = ipt_data2 + iptsize_data[0]; // rov
            E_Float *ipt_data4 = ipt_data3 + iptsize_data[0];
            E_Float *ipt_data5 = ipt_data4 + iptsize_data[0]; // roe
            E_Float *ipt_data6 = ipt_data5 + iptsize_data[0]; // ronutilde

            bvbs_inflow_supersonic_fich_(idir, lrhs_loc, nstep,
                                         neq_mtr, mobile_coef,
                                         param_int, ipt_ind_CL,
                                         param_real, c4, c5, c6,
                                         iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                                         ipt_data1, ipt_data2, ipt_data3, ipt_data4, ipt_data5, ipt_data6,
                                         iptsize_data[0], inc_bc);
          }
        }
        else if (bc_type == 3 || (bc_type == 4 && param_int[IFLOW] == 1) || bc_type == 5)
        {
          DEFAULT_STATE("BCWallInviscid")
          E_Float mobile_coef = 1.;
          if (nbdata != 0)
            mobile_coef = ipt_data[0];

          bvbs_wall_inviscid_(idir, lrhs_loc, nstep,
                              neq_mtr, mobile_coef,
                              param_int, ipt_ind_CL,
                              param_real, c4, c5, c6,
                              iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 31)
        {
          DEFAULT_STATE("BCWallExchange")
#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1;
          E_Int shift_tab = 0;
          if (nbdata != 0 && iptsize_data[0] == 1)
          {
            mobile_coef = ipt_data[0];
            shift_tab = 1;
          }

          /*
          E_Float* iptAvgField    = ipt_data + shift_tab;
          E_Float* iptparam_wmles = ipt_data + iptsize_data[shift_tab ];

          E_Int size_data = iptsize_data[1]-iptsize_data[0];
          if( shift_tab ==0) size_data = iptsize_data[0];

          E_Int Nechant   = int(iptparam_wmles[0]);
          E_Int sizetab   =  size_data/(Nechant*5);
          printf("echant %d %d %d %d %d %d \n", nbdata, sizetab,Nechant, iptsize_data[0], iptsize_data[1], iptsize_data[2] );
          */

          bvbs_wallexchange_(idir, lrhs_loc, nstep,
                             neq_mtr, mobile_coef,
                             param_int, ipt_ind_CL,
                             param_real, c4, c5, c6,
                             iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 30)
        {
          DEFAULT_STATE("BCWallModel")
          E_Float mobile_coef = 1.;
          if (nbdata != 0)
            mobile_coef = ipt_data[0];

          if (lrhs_loc == 0)
          {
            bvbs_wallmodel_(idir, lrhs_loc, nstep,
                            neq_mtr, mobile_coef,
                            param_int, ipt_ind_CL,
                            param_real, c4, c5, c6,
                            iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
          }
          else
          {
            bvbs_wallmodel_rhs_(idir, lrhs_loc, nstep,
                                neq_mtr, mobile_coef,
                                param_int, ipt_ind_CL,
                                param_real, c4, c5, c6,
                                iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
          }
        }
        else if ((bc_type == 6 || bc_type == 4) && param_int[IFLOW] > 1)
        {
          DEFAULT_STATE("BCWallViscous")
          E_Float mobile_coef = 1.;
          if (nbdata != 0)
            mobile_coef = ipt_data[0];

          bvbs_wall_viscous_adia_(idir, lrhs_loc, nstep,
                                  neq_mtr, mobile_coef,
                                  param_int, ipt_ind_CL,
                                  param_real, c4, c5, c6,
                                  iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }
        else if (bc_type == 12 && param_int[IFLOW] > 1)
        {
          DEFAULT_STATE("BCWallViscous")
#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1.;
          E_Float *random_bc = NULL;
          E_Int size_data = 1;
          if (nbdata != 0)
          {
            if (iptsize_data[0] == 1)
            {
              mobile_coef = ipt_data[0];
              random_bc = ipt_data + 1;
              size_data = iptsize_data[1];
            }
            else
            {
              random_bc = ipt_data;
              if (nbdata == 2)
              {
                E_Float *tmp = ipt_data + iptsize_data[0];
                size_data = iptsize_data[0];
                mobile_coef = tmp[0];
              }
            }
          }

          bvbs_wall_viscous_transition_(idir, lrhs_loc, nstep,
                                        neq_mtr, mobile_coef,
                                        param_int, ipt_ind_CL,
                                        param_real, c4, c5, c6,
                                        iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                                        random_bc, iptsize_data[0], inc_bc);
        }
        else if (bc_type == 32 && param_int[IFLOW] > 1)
        {
          DEFAULT_STATE("BCWallViscousIsothermal")
#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1.;
          E_Int size_data = 1;
          E_Float *ipt_data1 = NULL;
          if (nbdata != 0)
          {
            if (iptsize_data[0] == 1)
            {
              mobile_coef = ipt_data[0];
              ipt_data1 = ipt_data + 1;
              size_data = iptsize_data[1];
            }
            else
            {
              ipt_data1 = ipt_data;
              if (nbdata == 2)
              {
                E_Float *tmp = ipt_data + iptsize_data[0];
                mobile_coef = tmp[0];
                size_data = iptsize_data[0];
              }
            }
          }
          bvbs_wall_viscous_isothermal_(idir, lrhs_loc, nstep,
                                        neq_mtr, mobile_coef,
                                        param_int, ipt_ind_CL,
                                        param_real, c4, c5, c6,
                                        iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                                        ipt_data1, iptsize_data[0], inc_bc);
        }

        else if (bc_type == 11)
        {
          DEFAULT_STATE("BCPeriodic")
          E_Float mobile_coef = 1.;
          bvbs_periodique_(idir, lrhs_loc, nstep,
                           neq_mtr, mobile_coef,
                           param_int, ipt_ind_CL119,
                           param_real, c4, c5, c6,
                           iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc);
        }

        else if (bc_type == 16)
        {

          DEFAULT_STATE("BCOutpres")
#include "BC/INCREMENT_BC.h"
          E_Float mobile_coef = 1.;
          E_Int size_data = 1;
          E_Float *ipt_data1 = NULL;
          if (nbdata != 0)
          {
            if (iptsize_data[0] == 1)
            {
              mobile_coef = ipt_data[0];
              ipt_data1 = ipt_data + 1;
              size_data = iptsize_data[1];
            }
            else
            {
              ipt_data1 = ipt_data;
              if (nbdata == 2)
              {
                E_Float *tmp = ipt_data + iptsize_data[0];
                mobile_coef = tmp[0];
                size_data = iptsize_data[0];
              }
            }
          }

          if (nbdata > 2)
          {
            printf("ERROR BCOutpres: inappropriate number of BC data \n");
            exit(err);
          }

          bvbs_outpres_(idir, lrhs_loc, nstep,
                        neq_mtr, mobile_coef,
                        param_int, ipt_ind_CL,
                        param_real, c4, c5, c6,
                        iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                        ipt_data1, iptsize_data[0], inc_bc);
        }
        else if (bc_type == 17 && lrhs_loc == 0)
        {

          DEFAULT_STATE("BCinj1")
#include "BC/INCREMENT_BC.h"

          E_Float mobile_coef = 1.;

          E_Float *ipt_data1 = ipt_data;
          E_Float *ipt_data2 = ipt_data1 + iptsize_data[0];
          E_Float *ipt_data3 = ipt_data2 + iptsize_data[0];
          E_Float *ipt_data4 = ipt_data3 + iptsize_data[0];
          E_Float *ipt_data5 = ipt_data4 + iptsize_data[0];
          E_Float *ipt_data6 = ipt_data5 + iptsize_data[0];

          if (nbdata != 6 and param_int[IFLOW] == 3)
          {
            printf("ERROR BCInj1: inappropriate number of BC data \n");
            exit(err);
          }
          if (nbdata != 5 and param_int[IFLOW] <= 2)
          {
            printf("ERROR BCInj1: inappropriate number of BC data \n");
            exit(err);
          }

          bvbs_inflow_newton_(idir, lrhs_loc, nstep,
                              neq_mtr, mobile_coef,
                              param_int, ipt_ind_CL,
                              param_real, c4, c5, c6,
                              iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptmut, ipt_data_loc,
                              ipt_data1, ipt_data2, ipt_data3, ipt_data4, ipt_data5, ipt_data6,
                              iptsize_data[0], inc_bc);
        }

      } // if skip_loc

      if (lrhs == 2)
      {
        E_Int range = 2;
        E_Int without_ghost = 0;
        E_Float signe = -1.;
        // E_Int ipt_ind_CLgmres[6];

        for (E_Int i = 0; i < 6; i++)
          ipt_ind_CLgmres[i] = ipt_ind_CL[i];

        if (idir % 2 == 1)
        {
          ipt_ind_CLgmres[idir - 1] += range * without_ghost;
          ipt_ind_CLgmres[idir] += range;
        }
        else
        {
          range *= -1;
          ipt_ind_CLgmres[idir - 2] += range;
          ipt_ind_CLgmres[idir - 1] += range * without_ghost;
        }

        pre_bc_(param_int, signe, ipt_ind_CLgmres, iptrop, iptrop_gmres);
      }

    } // if skip_ipara
  } // fin boucle CL

  return err;
}
