     //
     //calcul norme L2 et Loo du RHS
     E_Int nbtask = ipt_omp[nitcfg-1]; 
     E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

     //printf("nbtaaask %d %d \n", nbtask, ptiter);
     E_Int barrier_residu = 0;
     for (E_Int ntask = 0; ntask < nbtask; ntask++)
       {
         E_Int pttask     = ptiter + ntask*(6+Nbre_thread_actif*7);
         E_Int nd         = ipt_omp[ pttask ];
         E_Int nd_subzone = ipt_omp[ pttask + 1 ];


         E_Int* ipt_ind_dm_loc = ipt_ind_dm[nd] + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
         E_Int* ipt_nidom_loc  = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
         E_Int  nb_subzone     = ipt_nidom_loc [nitcfg-1]; 

         shift_zone=0; shift_coe=0;
         for (E_Int n = 0; n < nd; n++)
           {
            shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
            shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
           }

          //
          //verrou rhs
          //
          E_Int type = 2;
          for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
          { 
               E_Int* verrou_lhs_thread= verrou_lhs + ntask*Nbre_thread_actif +th;
               verrou_c_( verrou_lhs_thread, type );
          }

          E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX];
          if(nitcfg==1)
          {
            size = param_int[nd][NEQ_COE]*param_int[nd][NDIMDX];
            flush_real_( size , iptcoe + shift_coe);
          }
          //sortie de la carte d residu du Newton
          //
          //
          //
          //Investiger pourquoi resultat depend de SIZE_SSDOM
#include  "Compute/residus_navier.h"

          if(lssiter_verif ==1  && nd_subzone ==0 && omp_mode==1 && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
          {
            E_Int type = 2;
           for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
             { 
               E_Int* verrou_lhs_thread= verrou_lhs + (mx_nidom + ntask)*Nbre_thread_actif +th;
               verrou_c_( verrou_lhs_thread, type );
             }
          }
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
#include "FastC/HPC_LAYER/OMP_MODE_BEGIN.h"
	    E_Float* increment = iptdrodm + shift_zone;
	    if (param_int[nd][NB_RELAX] > 1) increment = iptssortmp[nd];

            mjr_prim_from_cons_(param_int[nd], param_real[nd], param_real[nd]+VISCO, param_real[nd]+SA_REAL, ipt_ind_dm_thread, iptro_CL[nd], iptro_ssiter[nd], increment);
#include "FastC/HPC_LAYER/OMP_MODE_END.h"
      }
 }//if lexit_lu


