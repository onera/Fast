     //
     //calcul norme L2 et Loo du RHS
     shift_zone=0; shift_coe=0; nd_current=0;
     for (E_Int nd = 0; nd < nidom; nd++)
     {
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
          //
          //verrou rhs
          //
          E_Int type = 2;
          for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
          { 
               E_Int* verrou_lhs_thread= verrou_lhs + nd_current*Nbre_thread_actif +th;
               verrou_c_( verrou_lhs_thread, type );
          }

          E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX];
          if(nitcfg==1)
          {
            size = param_int[nd][NEQ_COE]*param_int[nd][NDIMDX];
            flush_real_( size , iptcoe + shift_coe);
          }
          //sortie de la carte d residu du Newton
#include  "Compute/residus_navier.h"
          nd_current +=1;

#include  "HPC_LAYER/OMP_MODE_END.h"
     shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
     shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
     }// fin loop zones


if(lexit_lu == 0 )
 {
    E_Int num_max_vect = param_int[0][NB_KRYLOV];
    E_Int nb_restart   = param_int[0][NB_RESTART];
    bool  continue_gmres(true);
    E_Int restart = 0;

    //Voir https://web.stanford.edu/class/cme324/saad.pdf pour les notations

    //VectG = Beta x e1, ayant subi les rotations de Givens
    //VectY est la sol de Hessenberg x VectY = VectG
    //Avec Hessenberg après rotation de Givens et on ignore la derniere ligne
    //de Hessengerg + VectG
    //Taille d'Hessenberg = (num_max_vect, num_max_vect - 1)
    //

    //while (continue_gmres && (restart < nb_restart))
    while (restart < nb_restart)
    {
#include  "Compute/Linear_solver/restart_gmres.h"
      restart++;
    }//loop restart

   //mise a jour nouvelle solution
    if (param_int[0][NB_RELAX] == 0)
      {
        shift_zone =0; nd_current =0;
        for (E_Int nd = 0; nd < nidom; nd++)
          {
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
	    E_Float* increment = iptdrodm + shift_zone;
	    if (param_int[nd][NB_RELAX] > 1) increment = iptssortmp[nd];

            mjr_prim_from_cons_(param_int[nd], param_real[nd], param_real[nd]+VISCO, param_real[nd]+SA_REAL, ipt_ind_dm_thread, iptro_CL[nd], iptro_ssiter[nd], increment);
            nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"

            shift_zone  = shift_zone  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
          }//loop zone
      }
 }//if lexit_lu


