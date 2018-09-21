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
#include "fastS.h"
#include "param_solver.h"
#if defined( _OPENMP )
#include <omp.h>
#endif

#define DEFAULT_STATE( string )                                                                   \
    FldArrayF data( 6 );                                                                          \
    if ( nbdata == 0 ) {                                                                          \
        ipt_data    = data.begin( );                                                              \
        ipt_data[0] = param_real[ROINF];                                                          \
        ipt_data[1] = param_real[VXINF] * param_real[ROINF];                                      \
        ipt_data[2] = param_real[VYINF] * param_real[ROINF];                                      \
        ipt_data[3] = param_real[VZINF] * param_real[ROINF];                                      \
        ipt_data[4] = param_real[PINF] / ( param_real[GAMMA] - 1. ) +                             \
                      0.5 * ( param_real[VXINF] * ipt_data[1] + param_real[VYINF] * ipt_data[2] + \
                              param_real[VZINF] * ipt_data[3] );                                  \
        ipt_data[5] = param_real[RONUTILDEINF];                                                   \
    }

using namespace std;
using namespace K_FLD;

//=============================================================================
E_Int K_FASTS::BCzone_d(
    E_Int& nd, E_Int& lrhs, E_Int& lcorner,
    // E_Int& ithread,  E_Int& Nbre_thread_actif, E_Int* ipt_thread_topology, E_Float* vteta, E_Float* roteta,
    E_Int* param_int, E_Float* param_real, E_Int& npass, E_Int* ipt_ind_dm, E_Int* ipt_ind_dm_thread, E_Int* ipt_ind_CL,
    E_Int* ipt_ind_CL119, E_Int* ipt_ind_CLgmres, E_Int* ishift_lu, E_Float* iptrop, E_Float* ipti, E_Float* iptj, E_Float* iptk, E_Float* iptx,
    E_Float* ipty, E_Float* iptz, E_Float* iptventi, E_Float* iptventj, E_Float* iptventk, E_Float* iptrop_d ) {
    E_Int err = 1;
    E_Int neq_mtr;
    E_Int itest[6], lskip[6], ind_rhs[6], ind_avg[6], ind_avg_thread[6], ind_mjr[6], ind_mjr_thread[6];

    E_Int lrhs_loc = 0;
    if (lrhs == 1) lrhs_loc = 1;

    itest[0] = 1;
    itest[1] = param_int[IJKV];
    itest[2] = 1;
    itest[3] = param_int[IJKV + 1];
    itest[4] = 1;
    itest[5] = param_int[IJKV + 2];

    for ( E_Int i = 0; i < 6; i++ ) lskip[i] = 1;

    E_Int* ipt_nijk = param_int + NIJK;
    E_Int* ipt_ijkv = param_int + IJKV;

    E_Float* iptijk;
    E_Float* ipventijk;

    E_Int pt_bcs = param_int[PT_BC];
    E_Int nb_bc  = param_int[ pt_bcs ];


    for ( E_Int ndf = 0; ndf < nb_bc; ndf++ ) {
        E_Int pt_bc = param_int[ pt_bcs + 1 + ndf];
        E_Int idir  = param_int[pt_bc + BC_IDIR];

        E_Int ipara = idir - 1;
        // on test si le sous domaine touche le bord du domaine et si la CL est bonne candidate a implicitation
        if ( lrhs_loc == 1 ) {
                          if ( ( ipt_ind_dm_thread[ipara] == itest[ipara] ) &&
                             ( ( param_int[pt_bc + BC_TYPE] >= 3 && param_int[pt_bc + BC_TYPE] <= 7 ) ||
                                 param_int[pt_bc + BC_TYPE] == 12 ) ) 
                           {
                             lskip[ipara] = 0;
                             lcorner      = 0;
                           }
	}
        else {lskip[ipara] = 1;}

    }  // fin loop CL warmup

    if ( lrhs_loc == 1 ) {
        E_Int ificmax = 1;
        ishift_lu[0]  = ipt_ind_dm_thread[0] - ( 1 - lskip[0] ) * ificmax;
        ishift_lu[1]  = ipt_ind_dm_thread[1] + ( 1 - lskip[1] ) * ificmax;
        ishift_lu[2]  = ipt_ind_dm_thread[2] - ( 1 - lskip[2] ) * ificmax;
        ishift_lu[3]  = ipt_ind_dm_thread[3] + ( 1 - lskip[3] ) * ificmax;
        ishift_lu[4]  = ipt_ind_dm_thread[4] - ( 1 - lskip[4] ) * ificmax;
        ishift_lu[5]  = ipt_ind_dm_thread[5] + ( 1 - lskip[5] ) * ificmax;

        E_Int idirmax                           = 6;
        if ( param_int[ITYPZONE] == 3 ) idirmax = 4;

        // extrapolation du RHS par defaut sur toutes les faces ou on etend la resolution LU.  (Sinon probleme si face
        // domaine= paroi + raccord)
        for ( E_Int idir = 0; idir < idirmax; idir++ ) {
            for ( E_Int l = 0; l < 6; l++ ) ind_rhs[l] = ipt_ind_dm_thread[l];
            E_Int       lskip_loc;
            if ( idir == 0 ) {
                ind_rhs[0] = 1;
                ind_rhs[1] = 1;
            }
            if ( idir == 2 ) {
                ind_rhs[2] = 1;
                ind_rhs[3] = 1;
            }
            if ( idir == 4 ) {
                ind_rhs[4] = 1;
                ind_rhs[5] = 1;
            }
            if ( idir == 1 ) { ind_rhs[0] = ind_rhs[1]; }
            if ( idir == 3 ) { ind_rhs[2] = ind_rhs[3]; }
            if ( idir == 5 ) { ind_rhs[4] = ind_rhs[5]; }

            if ( lskip[idir] == 0 ) {
                E_Int idir_loc   = idir + 1;
                E_Int typebc_loc = 0;
                indice_cl_sdm_( idir_loc, npass, lskip_loc, typebc_loc, ipt_nijk[3], ipt_nijk[4],  // IN
                                ipt_ijkv, ind_rhs, ipt_ind_dm, ipt_ind_dm_thread,                  // IN
                                ipt_ind_CL, ipt_ind_CL119 );

                E_Int eq_deb = 1;
                if ( lskip_loc == 0 )
                    bvbs_extrapolate_( idir_loc, lrhs_loc, eq_deb, param_int, ipt_ind_CL119, param_real[RONUTILDEINF],
                                       iptrop );
            }
        }
    }
    /*
      for (E_Int ndf = 0; ndf < nb_bc; ndf++)
      {
        E_Int pt_bc  = param_int[BC_NBBC + 1 + ndf];
        E_Int idir   = param_int[ pt_bc + BC_IDIR];
        E_Int nbdata = param_int[ pt_bc + BC_NBDATA];

        E_Int* iptsize_data = param_int +  pt_bc + BC_NBDATA +1;


        E_Int bc_type       = param_int[pt_bc + BC_TYPE];

        E_Float* ipt_data;
        if (nbdata != 0) ipt_data = param_real + param_int[BC_NBBC + 1 + ndf + nb_bc];

        E_Int ipara  = idir -1;
        if (bc_type == 16 && lrhs==0)
            {
               E_Int flag_avg = 2; // moyenne azymutal suivant j
               E_Int* ind_fen =  param_int+pt_bc + BC_FEN;
               ind_mjr[0:6]   = ind_fen[0:6];
               ind_avg[0:6]   = ind_fen[0:6];

               if (idir <=2)
                  { if      (flag_avg==2){ind_avg[2] = ind_fen[2]; ind_avg[3]= ind_fen[2]; ind_mjr[4] = ind_fen[4];
    ind_mjr[5]= ind_fen[4];}
                    else if (flag_avg==3){ind_avg[4] = ind_fen[4]; ind_avg[5]= ind_fen[4]; ind_mjr[2] = ind_fen[2];
    ind_mjr[3]= ind_fen[2];}
                  }

               E_Int  inc_bc;
               if  ( idir <=2) { inc_bc = ind_fen[3] - ind_fen[2] +1;}  // nombre element de la fenetre dans la
    direction J
               else            { inc_bc = ind_fen[1] - ind_fen[0] +1;}  // nombre element de la fenetre dans la
    direction I

               E_Int itypcp_loc = 2;
               indice_boucle_lu_(nd, ithread, Nbre_thread_actif, itypcp_loc,
                               ind_avg,
                               ipt_thread_topology, ind_avg_thread);
               indice_boucle_lu_(nd, ithread, Nbre_thread_actif, itypcp_loc,
                               ind_mjr,
                               ipt_thread_topology, ind_mjr_thread);

               if    (  flag_avg==2) {ind_avg_thread[2] = ind_fen[2]; ind_avg_thread[3]= ind_fen[3]; ind_mjr_thread[4] =
    ind_fen[4]; ind_mjr_thread[5]= ind_fen[5];}

           bvbs_updatepressure_(idir, ithread, ind_avg_thread, ind_mjr_thread, param_int, ipt_ind_CL, param_real,
                                    ipty, iptz, iptrop, ipt_data, iptsize_data[0], vteta, roteta, inc_bc );
            }
      }
    #pragma omp barrier
    */

    for ( E_Int ndf = 0; ndf < nb_bc; ndf++ ) {
        E_Int pt_bc  = param_int[pt_bcs+ 1 + ndf];

        E_Int idir   = param_int[pt_bc + BC_IDIR];
        E_Int nbdata = param_int[pt_bc + BC_NBDATA];
        E_Int bc_type= param_int[pt_bc + BC_TYPE];

        //printf("ptbc= %d , idir= %d , nbdat= %d , bctype=  %d \n", pt_bc, idir,nbdata,  bc_type );

        E_Int* iptsize_data = param_int + pt_bc + BC_NBDATA + 1;

        E_Int ipara = idir - 1;

        if ( lskip[ipara] == 0 || lrhs_loc == 0 )  // test pour savoir si face ipara est skipper en LU pour le rhs
        {
            if ( idir <= 2 ) {
                iptijk    = ipti;
                neq_mtr   = param_int[NEQ_IJ];
                ipventijk = iptventi;
            } else if ( idir <= 4 ) {
                iptijk    = iptj;
                neq_mtr   = param_int[NEQ_IJ];
                ipventijk = iptventj;
            } else {
                iptijk    = iptk;
                neq_mtr   = param_int[NEQ_K];
                ipventijk = iptventk;
            }

            // calcul intersection sous domaine et fenetre CL
            E_Int lskip_loc;
            indice_cl_sdm_( idir, npass, lskip_loc, param_int[pt_bc + BC_TYPE], ipt_nijk[3], ipt_nijk[4],  // IN
                            ipt_ijkv, param_int + pt_bc + BC_FEN, ipt_ind_dm, ipt_ind_dm_thread,           // IN
                            ipt_ind_CL, ipt_ind_CL119 );                                                   // OUT

            E_Int bc_type = param_int[pt_bc + BC_TYPE];
            E_Float* ipt_data;
	    
            if ( nbdata != 0 ) ipt_data = param_real + param_int[pt_bcs + 1 + ndf + nb_bc];

            if ( lskip_loc == 0 ) {
                // RANS:LES extrapolation
                E_Int eq_deb                             = 1;
                if ( bc_type == 14 && lrhs_loc == 0 ) eq_deb = 6;

                if ( lrhs_loc == 1 && ( bc_type == 0 || bc_type == 1 || bc_type == 2 || bc_type == 10 || bc_type == 11 ||
                                    bc_type == 14 || bc_type == 16 || bc_type == 17 ) )
                    bc_type = 0;  // si rhs, on extrapole sauf si wall

                if ( bc_type == 0 ) { bvbs_extrapolate_d_( idir, lrhs_loc, eq_deb, param_int, ipt_ind_CL119, param_real[RONUTILDEINF], iptrop, iptrop_d );
                
                } else if ( bc_type == 1) {
                    DEFAULT_STATE( "BCFarfield" )
                    // MUSCL
                    E_Float c4, c5, c6;
                    c4 = 5. / 6.;
                    c5 = 2. / 6.;
                    c6 = -1. / 6.;

                    bvbs_farfield_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6, ipventijk, iptijk, iptrop, iptrop_d, ipt_data );

                } else if ( bc_type == 10 ) {
                    DEFAULT_STATE( "BCOutflow" )
                    // MUSCL
                    E_Float c4, c5, c6;
                    c4 = 5. / 6.;
                    c5 = 2. / 6.;
                    c6 = -1. / 6.;

                    bvbs_outflow_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6, ipventijk, iptijk, iptrop, iptrop_d, ipt_data );
                } else if ( bc_type == 13 ) {
                    DEFAULT_STATE( "BCInflow" )
                    // MUSCL
                    E_Float c4, c5, c6;
                    c4 = 5. / 6.;
                    c5 = 2. / 6.;
                    c6 = -1. / 6.;

                    if (iptsize_data[0] <=6)
                     {
                      bvbs_inflow_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6, ipventijk, iptijk, iptrop, iptrop_d, ipt_data );
                     }
                    else
                     { 
#                      include "BC/INCREMENT_BC.h"

                       E_Int size_work = ipt_ind_CL[1] - ipt_ind_CL[0] + 1;

                       E_Float* ipt_data1 = ipt_data;
                       E_Float* ipt_data2 = ipt_data1 + iptsize_data[0];
                       E_Float* ipt_data3 = ipt_data2 + iptsize_data[0];
                       E_Float* ipt_data4 = ipt_data3 + iptsize_data[0];
                       E_Float* ipt_data5 = ipt_data4 + iptsize_data[0];
                       E_Float* ipt_data6 = ipt_data5 + iptsize_data[0];

                       bvbs_inflow_fich_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6, ipventijk,
                                            iptijk, iptrop, iptrop_d, ipt_data1, ipt_data2, ipt_data3, ipt_data4, ipt_data5,
                                            ipt_data6, iptsize_data[0], inc_bc, size_work );
                     }

                } else if ( bc_type == 2 ) {
                    DEFAULT_STATE( "BCSupersonic" )
                    // MUSCL
                    E_Float c4, c5, c6;
                    c4 = 5. / 6.;
                    c5 = 2. / 6.;
                    c6 = -1. / 6.;

                    bvbs_inflow_supersonic_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6,
                                               ipventijk, iptijk, iptrop, iptrop_d, ipt_data );
                } else if ( bc_type == 3 || ( bc_type == 4 && param_int[IFLOW] == 1 ) || bc_type == 5 ) {
                    E_Float mobile_coef            = 1.;
                    if ( nbdata != 0 ) mobile_coef = ipt_data[0];

                    bvbs_wall_inviscid_d_( idir, lrhs_loc, neq_mtr, mobile_coef, param_int, ipt_ind_CL, ipventijk, iptijk, iptrop, iptrop_d );

                } else if ( ( bc_type == 6 || bc_type == 4 ) && param_int[IFLOW] > 1 ) {
                    E_Float mobile_coef            = 1.;
                    if ( nbdata != 0 ) mobile_coef = ipt_data[0];

                    bvbs_wall_viscous_adia_d_( idir, lrhs_loc, neq_mtr, mobile_coef, param_int, ipt_ind_CL, ipventijk, iptijk,
                                               iptrop, iptrop_d );
                } else if ( bc_type == 12 && param_int[IFLOW] > 1 ) {
                    E_Float mobile_coef            = 1.;
                    if ( nbdata != 0 ) mobile_coef = ipt_data[0];

                    bvbs_wall_viscous_transition_d_( idir, lrhs_loc, neq_mtr, mobile_coef, param_int, ipt_ind_CL, param_real,
                                                     iptx, ipty, iptz, ipventijk, iptijk, iptrop, iptrop_d );
                }

                else if ( bc_type == 11 ) {
                    bvbs_periodique_d_( idir, lrhs_loc, param_int, ipt_ind_CL119, iptrop, iptrop_d );
                }
                /*else if (bc_type == 15 )
                    {
                   //if (ipt_data[0] != 0. or  ipt_data[1] != 0. or  ipt_data[2] != 0.)
                    // periodicite par rotation
                    // bvbs_periodique_azimuthal_(idir, lrhs, param_int, ipt_ind_CL, iptrop, ipt_data);
                   //else
                     // Voir Avec JCB la raison de cet appel
                     // parametres de rotation nuls => periodicite par translation
                     bvbs_periodique_(idir, lrhs, param_int, ipt_ind_CL119, iptrop);
                    }*/
                else if ( bc_type == 16 ) {
                    // MUSCL
                    E_Float c4, c5, c6;
                    c4 = 5. / 6.;
                    c5 = 2. / 6.;
                    c6 = -1. / 6.;

#                   include "BC/INCREMENT_BC.h"

                    bvbs_outpres_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6, ipventijk,
                                     iptijk, iptrop, iptrop_d, ipt_data, iptsize_data[0], inc_bc );

                } else if ( bc_type == 17 && lrhs_loc == 0 ) {
                    // MUSCL
                    E_Float c4, c5, c6;
                    c4 = 5. / 6.;
                    c5 = 2. / 6.;
                    c6 = -1. / 6.;

#                   include "BC/INCREMENT_BC.h"

                    E_Int size_work = ipt_ind_CL[1] - ipt_ind_CL[0] + 1;

                    E_Float* ipt_data1 = ipt_data;
                    E_Float* ipt_data2 = ipt_data1 + iptsize_data[0];
                    E_Float* ipt_data3 = ipt_data2 + iptsize_data[0];
                    E_Float* ipt_data4 = ipt_data3 + iptsize_data[0];
                    E_Float* ipt_data5 = ipt_data4 + iptsize_data[0];
                    E_Float* ipt_data6 = ipt_data5 + iptsize_data[0];

                    bvbs_inflow_newton_d_( idir, lrhs_loc, neq_mtr, param_int, ipt_ind_CL, param_real, c4, c5, c6, ipventijk,
                                           iptijk, iptrop, iptrop_d, ipt_data1, ipt_data2, ipt_data3, ipt_data4, ipt_data5,
                                           ipt_data6, iptsize_data[0], inc_bc, size_work );

                } 
            }  // if skip_loc
        }      // if skip_ipara
    }          // fin boucle CL

    return err;
}
