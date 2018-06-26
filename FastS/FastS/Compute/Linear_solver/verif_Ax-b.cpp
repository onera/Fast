E_Int sizevect = 0;
for (E_Int nd = 0; nd < nidom; ++nd)
    sizevect += param_int[nd][NEQ] * param_int[nd][NDIMDX];

E_Float* ipt_tmpverif = new E_Float[sizevect]();

// shift_zone = 0;
// for (E_Int nd = 0; nd < nidom; nd++)
//   {
//     E_Float* krylov_in    = ipt_tmpverif + shift_zone;
//     for (E_Int i = 0; i < param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ]; i++)
//       krylov_in[i] = iptdrodm[i];

//     shift_zone  = shift_zone  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
//   }

shift_zone = 0;
for (E_Int nd = 0; nd < nidom; nd++)
  {
    //A*V_kr-1 tapenade soon   //sur ind_sdm
    E_Float* krylov_in = iptdrodm + shift_zone;
    E_Float* krylov_out= ipt_tmpverif + shift_zone;
    iptkrylov_transfer[nd] = krylov_out;
#include "HPC_LAYER/OMP_MODE_BEGIN.h"

    dp_dw_vect_(param_int[nd], param_real[nd],ipt_ind_dm_thread , iptro_ssiter[nd], krylov_in,  krylov_out);

#include "HPC_LAYER/OMP_MODE_END.h"
    shift_zone  = shift_zone  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
  }//loop zone

#pragma omp master
    { //Raccord V0
      setInterpTransfersFastS(iptkrylov_transfer, ndimdx_transfer, param_int_tc,
    			      param_real_tc, param_int_ibc, param_real_ibc, param_real[0][PRANDT],
         		      it_target, nidom, ipt_timecount, mpi);
    }
#pragma omp barrier
    E_Int lrhs = 0; E_Int lcorner = 0; E_Int npass = 0; E_Int ipt_shift_lu[6];
    nd_current=0;shift_zone =0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
       E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
       E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
       E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
       E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;
#include  "HPC_LAYER/OMP_MODE_BEGIN.h"

       E_Int ierr = BCzone_d(nd, lrhs , lcorner, param_int[nd], param_real[nd], npass,
			     ipt_ind_dm_loc, ipt_ind_dm_thread,
			     ipt_ind_CL_thread, ipt_ind_CL119 ,  ipt_ind_CLgmres, ipt_shift_lu,
			     iptro_ssiter[nd],
			     ipti[nd]     , iptj[nd]            , iptk[nd],
			     iptx[nd]     , ipty[nd]            , iptz[nd],
			     iptventi[nd] , iptventj[nd]        , iptventk[nd],
			     iptkrylov_transfer[nd]);

       //correct_coins_(nd,  param_int[nd], ipt_ind_dm_thread , iptkrylov_transfer[nd]);

            nd_current +=1;
#include    "HPC_LAYER/OMP_MODE_END.h"
	    shift_zone  = shift_zone  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
      }//loop zone

#pragma omp barrier

    shift_zone=0;
    shift_wig =0;
    shift_coe =0;
    nd_current=0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
       E_Float* krylov_in    = iptkrylov[nd] +  3    * param_int[nd][NEQ] * param_int[nd][NDIMDX];
       E_Float* rop_ssiter_d = iptkrylov_transfer[nd];

       E_Int lmin = 10;
       if (param_int[nd][ITYPCP] == 2) lmin = 4;

#include "Compute/Linear_solver/dRdp_tapenade_test.cpp"

	shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
	shift_wig  = shift_wig  + param_int[nd][ NDIMDX ]*3;
	shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
      }

    shift_zone=0;
    nd_current=0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
        E_Float* krylov_in = iptdrodm + shift_zone;
        E_Float* krylov_out= ipt_tmpverif + shift_zone;
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
             id_vect_(param_int[nd], ipt_ind_dm_thread, ipt_drodmd + shift_zone, krylov_out, krylov_in);
             nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"

        shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
      }

cout << "Valeur norme b = " << save << endl;

normL2 = 0., normL2_sum = 0.;
nd_current =0; shift_zone = 0;
for (E_Int nd = 0; nd < nidom; nd++)
  {
    E_Float* krylov_i   = iptkrylov[nd];
    E_Float* krylov_krp1= ipt_tmpverif + shift_zone;

#include "HPC_LAYER/OMP_MODE_BEGIN.h"
    scal_prod_(param_int[nd], ipt_ind_dm_thread, krylov_krp1, krylov_krp1, normL2);

    vect_rvect_(param_int[nd], ipt_ind_dm_thread, krylov_krp1, krylov_i, save, normL2_sum);

    nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
    shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
  }

cout << "Valeur norme Ax = " << sqrt(normL2) << endl;
cout << "Valeur norme Ax - b = " << sqrt(normL2_sum) << endl;

delete[] ipt_tmpverif;
