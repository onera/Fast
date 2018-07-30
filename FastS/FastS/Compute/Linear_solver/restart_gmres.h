  E_Float value = 0.; E_Float sum_value = 0.; E_Float normL2 = 0., normL2_sum = 0.; E_Int mjrnewton = 0;

  E_Float* ipt_ssor_shift; E_Int indice;
  // norm L2 krylov pour openmp
  for (E_Int th = 0; th < Nbre_thread_actif; th++) { normL2_sum += ipt_norm_kry[th];}

  //Normalisation de V0 + initialisation de G = Beta x e1
  normL2_sum = sqrt(normL2_sum);

#pragma omp single
  {
     ipt_VectG[0] = normL2_sum;
  }

  E_Float epsi_linear = normL2_sum*param_real[0][EPSI_LINEAR];
//E_Float epsi_linear = param_real[0][EPSI_LINEAR];
  //if (ithread==1) printf(" norm Vo %g \n", normL2_sum);

  // Pour verif Ax-b
  E_Float save = normL2_sum;

E_Float lhs_beg = omp_get_wtime();


    //iptkrylov[nd] = iptkrylov[nd] / normL2_sum
    nd_current =0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
#include "HPC_LAYER/OMP_MODE_BEGIN.h"

	  normalisation_vect_(normL2_sum, param_int[nd], ipt_ind_dm_thread , iptkrylov[nd]);

          nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
      }
E_Float lhs_end = omp_get_wtime();
if (ithread==1) printf("cpu normalisation_vect %g \n", lhs_end-lhs_beg);
//
//
// Loop sur vectuer krylov
//
//
E_Int kr = 0;
while ((kr < num_max_vect - 1) && continue_gmres)
  {
#pragma omp barrier
     //if( ithread==1) printf("iter krylov =  %d  \n" , kr);

lhs_beg = omp_get_wtime();
    // 2.1) Calcul de V_kr = A * V_kr-1
        shift_coe=0; shift_zone=0; nd_current=0;
	for (E_Int nd = 0; nd < nidom; nd++)
	  {
	   //A*V_kr-1 tapenade soon   //sur ind_sdm
	   E_Float* krylov_in = iptkrylov[nd] +  kr     *param_int[nd][NEQ] * param_int[nd][NDIMDX];
	   E_Float* krylov_out= iptkrylov[nd] + (kr + 1)*param_int[nd][NEQ] * param_int[nd][NDIMDX];

	   //E_Float*       ssor=  ipt_ssor     +     shift_zone;

#include "HPC_LAYER/OMP_MODE_BEGIN.h"

#include "Compute/LU/prep_lussor.h"
	   indice = nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + ithread - 1;

	   invlu_(nd                     , nitcfg      ,nitrun, param_int[nd], param_real[nd],
	   	  ipt_ind_dm_thread      , ipt_ind_dm_thread       , mjrnewton             ,
	   	  iptrotmp[nd]           , iptro_ssiter[nd]        , krylov_in             , ipt_gmrestmp[nd],
	   	  ipti[nd]               , iptj[nd]                , iptk[nd]              ,
	   	  iptventi[nd]           , iptventj[nd]            , iptventk[nd]          ,
	   	  iptcoe  + shift_coe    , ipt_ssor_shift          , iptssortmp[nd]        , ipt_ssor_size[ indice ]);

	   if (param_int[nd][NB_RELAX] == 1) krylov_in = ipt_gmrestmp[nd];
	   else if (param_int[nd][NB_RELAX] > 1) krylov_in = iptssortmp[nd];
	   
             dp_dw_vect_(param_int[nd], param_real[nd],ipt_ind_dm_thread , iptro_ssiter[nd], krylov_in,  krylov_out);

             nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
	   shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
	   shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ     ];

	   //Tableau de pointeur pour les raccords V_kr
	   iptkrylov_transfer[nd] = krylov_out;

	  }//loop zone
lhs_end = omp_get_wtime();
if (ithread==1) printf("cpu invlu + dpdw       %g \n", lhs_end-lhs_beg);

#pragma omp barrier
    //Reinitialisation verrou omp
    //
    nd_current=0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
#include  "HPC_LAYER/OMP_MODE_BEGIN.h"
           E_Int l =  nd_current*mx_synchro*Nbre_thread_actif  + (ithread_loc-1)*mx_synchro;
           for (E_Int i = 0;  i < mx_synchro ; i++) { ipt_lok[ l + i ]  = 0; }
           nd_current +=1;
#include  "HPC_LAYER/OMP_MODE_END.h"
      }//loop zone


    //
    //
    // fillghostCell
    //
    //
lhs_beg = omp_get_wtime();
#pragma omp master
    { //Raccord V0
      setInterpTransfersFastS(iptkrylov_transfer, ndimdx_transfer, param_int_tc,
			      param_real_tc, param_int_ibc, param_real_ibc, param_real[0][PRANDT],
         		      it_target, nidom, ipt_timecount, mpi);
    }

lhs_beg = omp_get_wtime();
#pragma omp barrier
lhs_end = omp_get_wtime();
if (ithread==1) printf("transfert              %g \n", lhs_end-lhs_beg);

    E_Int lrhs = 2; E_Int lcorner = 0; E_Int npass = 0; E_Int ipt_shift_lu[6];
    nd_current=0;
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
			         iptkrylov_transfer[nd] );

/*            E_Int ierr = BCzone(nd, lrhs , lcorner, param_int[nd], param_real[nd], npass,
	           		 ipt_ind_dm_loc, ipt_ind_dm_thread, 
			         ipt_ind_CL_thread, ipt_ind_CL119 ,  ipt_ind_CLgmres, ipt_shift_lu,
			         iptro_ssiter[nd],
		    	         ipti[nd]     , iptj[nd]            , iptk[nd],
			         iptx[nd]     , ipty[nd]            , iptz[nd],
			         iptventi[nd] , iptventj[nd]        , iptventk[nd],
			         iptkrylov_transfer[nd]);
*/

            correct_coins_(nd,  param_int[nd], ipt_ind_dm_thread , iptkrylov_transfer[nd]);

            nd_current +=1;
#include    "HPC_LAYER/OMP_MODE_END.h"
      }//loop zone
lhs_end = omp_get_wtime();
if (ithread==1) printf("bczone                 %g \n", lhs_end-lhs_beg);


lhs_beg = omp_get_wtime();
#pragma omp barrier

    //
    //
    // produit matrice:vecteur.
    // In: rop_ssiter, rop_ssiter_d
    // Out: drodmd
    shift_zone    =0;
    shift_wig     =0;
    shift_coe     =0;
    E_Int shift_mu=0;
    nd_current    =0;
    E_Int mjr_dt  =0;
    ipt_norm_kry[ithread-1]=0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
       E_Float* krylov_in    = iptkrylov[nd] +  kr    * param_int[nd][NEQ] * param_int[nd][NDIMDX];
       E_Float* rop_ssiter_d = iptkrylov[nd] + (kr+1) * param_int[nd][NEQ] * param_int[nd][NDIMDX];

       E_Int lmin = 10;
       if (param_int[nd][ITYPCP] == 2) lmin = 4;
       
     
#include "Compute/Linear_solver/dRdp_tapenade.cpp"

	shift_mu   = shift_mu   + param_int[nd][ NDIMDX ];
	shift_wig  = shift_wig  + param_int[nd][ NDIMDX ]*3;
	shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
	shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
      }
#pragma omp barrier
lhs_end = omp_get_wtime();
if (ithread==1) printf("tapenade               %g \n", lhs_end-lhs_beg);

/*
    // norm L2 pour verrif bug
    for (E_Int th = 0; th < Nbre_thread_actif; th++) { normL2_sum += ipt_norm_kry[th];}
    normL2_sum = sqrt(normL2_sum);
    if(ithread==1) printf("tapenade   %f %d  \n",normL2_sum , ithread);
*/

lhs_beg = omp_get_wtime();
    // kry(kr+1) = kry(kr) + drodmd
    shift_zone=0;
    nd_current=0;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
        E_Float* krylov_in = iptkrylov[nd] +  kr    * param_int[nd][NEQ] * param_int[nd][NDIMDX];
        E_Float* krylov_out= iptkrylov[nd] + (kr+1) * param_int[nd][NEQ] * param_int[nd][NDIMDX];
	if (param_int[nd][NB_RELAX] == 1) krylov_in = ipt_gmrestmp[nd];
	else if (param_int[nd][NB_RELAX] > 1) krylov_in = iptssortmp[nd];
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
             id_vect_(param_int[nd], ipt_ind_dm_thread, ipt_drodmd + shift_zone, krylov_out, krylov_in);
             nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"

        shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
      }
lhs_end = omp_get_wtime();
if (ithread==1) printf("id_vect                %g \n", lhs_end-lhs_beg);

lhs_beg = omp_get_wtime();

    // 2.2) Orthonormalisation du vecteur V_kr
    //
    //
    for (E_Int i = 0; i < kr + 1; i++)
      {
//lhs_beg = omp_get_wtime();
        // Produit scalaire consecutif de Gram Schmidt modifie
        // de V_0 à V_kr-1
        // value =  V_ki . V_kr-1   sur ind-sdm (pas de dependance i+1)
        //
	sum_value              = 0.;
	normL2_sum             = 0.;
        ipt_norm_kry[ithread-1]= 0.;
        nd_current             = 0;
	for (E_Int nd = 0; nd < nidom; nd++)
	  {
	   E_Float* krylov_i   = iptkrylov[nd] +   i    * param_int[nd][NEQ] * param_int[nd][NDIMDX];
	   E_Float* krylov_krp1= iptkrylov[nd] + (kr+1) * param_int[nd][NEQ] * param_int[nd][NDIMDX];

#include "HPC_LAYER/OMP_MODE_BEGIN.h"

	       scal_prod_(param_int[nd], ipt_ind_dm_thread, krylov_krp1, krylov_i, ipt_norm_kry[ithread-1]);

               nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
	  }//loop zone
//lhs_end = omp_get_wtime();
//if (ithread==1) printf("scal_prod              %g \n", lhs_end-lhs_beg);
//lhs_beg = omp_get_wtime();

#pragma omp barrier
        // norm L2 krylov pour openmp
        for (E_Int th = 0; th < Nbre_thread_actif; th++) { sum_value += ipt_norm_kry[th];}

     //if(ithread==1)printf("orthonorm  %g %d  \n",sum_value , ithread);

        //Affectation des produits scalaires dans la matrice d'Hessenberg
#pragma omp single
        {
  	 E_Float* Hessenberg_i = ipt_Hessenberg + i * (num_max_vect - 1);
	 Hessenberg_i[kr]= sum_value;
        }

        //A chaque produit scalaire (V_kr, V_i) on fait
        //V_kr = V_kr - (V_kr, V_i)V_i
        //V_kr = V_kr - sum_value*V_i   sur ind-sdm (pas de dependance i+1)
        //
        // + calcul de la norme L2^2 mais seulement la derniere est utilisee
        ipt_norm_kry[ithread-1]= 0.;
        nd_current =0;
	for (E_Int nd = 0; nd < nidom; nd++)
	  {
	   E_Float* krylov_i   = iptkrylov[nd] +   i    * param_int[nd][NEQ] * param_int[nd][NDIMDX];
	   E_Float* krylov_krp1= iptkrylov[nd] + (kr+1) * param_int[nd][NEQ] * param_int[nd][NDIMDX];

#include "HPC_LAYER/OMP_MODE_BEGIN.h"
	       vect_rvect_(param_int[nd], ipt_ind_dm_thread, krylov_krp1, krylov_i, sum_value, ipt_norm_kry[ithread-1]);

               nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
          }
//lhs_end = omp_get_wtime();
//if (ithread==1) printf("vect_rvect             %g %d \n", lhs_end-lhs_beg, i);

      } // loop i< kr
        //
        //
        //
        //
lhs_end = omp_get_wtime();
if (ithread==1) printf("vect_rvect             %g %d \n", lhs_end-lhs_beg, kr);


#pragma omp barrier

    // norm L2 krylov pour openmp
    for (E_Int th = 0; th < Nbre_thread_actif; th++) { normL2_sum += ipt_norm_kry[th];}

    normL2_sum = sqrt(normL2_sum);
    //if (ithread==1) printf("norm rvect  %f %d  \n",normL2_sum , ithread);

    //Normalisation de V_kr + affectation de la norme sur la sous diagonale d'Hessenberg
#pragma omp single
    {
      E_Float tmp;
      E_Int shiftline = num_max_vect - 1;

      ipt_Hessenberg[ (kr+1)*shiftline + kr ] = normL2_sum;
      
      /* for (E_Int i = 0; i < kr + 2 ; i++) */
      /* 	{ */
      /* 	  E_Float* Hessenberg_i   = ipt_Hessenberg + i * (num_max_vect - 1); */
      /* 	  printf("hess AVANT \n " ); */
      /* 	  for (E_Int j = 0; j < kr + 1; j++) { printf(" %f %d %d ", Hessenberg_i[j], i,j ); } */
      /* 	  printf("\n" ); */
      /* 	} */

      for (E_Int i = 0; i < kr; i++)
      	{
          E_Int Iline   =  i    *shiftline;
          E_Int Ip1line = (i+1) *shiftline;

      	  tmp =  ipt_Hessenberg[ Iline + kr ];

          ipt_Hessenberg[ Iline   + kr ] =  ipt_givens[ i            ] * tmp + ipt_givens[ i+ shiftline ] * ipt_Hessenberg[ Ip1line + kr ];
          ipt_Hessenberg[ Ip1line + kr ] = -ipt_givens[ i+ shiftline ] * tmp + ipt_givens[ i            ] * ipt_Hessenberg[ Ip1line + kr ];
      	}

      tmp = sqrt( ipt_Hessenberg[ (kr+1)*shiftline + kr ]*ipt_Hessenberg[ (kr+1)*shiftline + kr ] + ipt_Hessenberg[ kr*(shiftline+1) ]*ipt_Hessenberg[ kr*(shiftline+1) ]);
      ipt_givens[ kr             ]  = ipt_Hessenberg[   kr  *(shiftline+1)     ]/ tmp;
      ipt_givens[ kr + shiftline ]  = ipt_Hessenberg[ (kr+1)* shiftline   + kr ]/ tmp;

      tmp = ipt_Hessenberg[ kr*(shiftline+1) ];
      ipt_Hessenberg[   kr  *(shiftline+1)   ] =   ipt_givens[ kr             ] * tmp + ipt_givens[ kr + shiftline ] *ipt_Hessenberg[ (kr+1)*shiftline + kr ];
      ipt_Hessenberg[ (kr+1)* shiftline + kr ] = - ipt_givens[ kr + shiftline ] * tmp + ipt_givens[ kr             ] *ipt_Hessenberg[ (kr+1)*shiftline + kr ];

      ipt_VectG[ kr + 1 ] = - ipt_givens[ kr + shiftline ] * ipt_VectG[ kr ];
      ipt_VectG[ kr ]     =   ipt_givens[ kr             ] * ipt_VectG[ kr ];

      /* for (E_Int i = 0; i < kr + 2 ; i++) */
      /* 	{ */
      /* 	  E_Float* Hessenberg_i   = ipt_Hessenberg + i * (num_max_vect - 1); */
      /* 	  printf("hess APRES \n " ); */
      /* 	  for (E_Int j = 0; j < kr + 1; j++) { printf(" %f %d %d ", Hessenberg_i[j], i,j ); } */
      /* 	  printf("\n" ); */
      /* 	} */

      //cout << "Residu GMRES = " << K_FUNC::E_abs(ipt_VectG[kr + 1]) << endl;
    }// end single
lhs_end = omp_get_wtime();
if (ithread==1) printf("Normalisation de V_kr  %g \n", lhs_end-lhs_beg);
lhs_beg = omp_get_wtime();

    nd_current =0;
    ipt_norm_kry[ithread-1]= 0.;
    for (E_Int nd = 0; nd < nidom; nd++)
      {
       E_Float* krylov_krp1= iptkrylov[nd] + (kr+1) * param_int[nd][NEQ] * param_int[nd][NDIMDX];
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
	       normalisation_vect_(normL2_sum, param_int[nd], ipt_ind_dm_thread , krylov_krp1);

             /*
               for ( E_Int i=0 ; i< kr+1; i++) 
               {
                  E_Float* krylov_i   = iptkrylov[nd] + i* param_int[nd][NEQ] * param_int[nd][NDIMDX];

                  ipt_norm_kry[ithread_loc-1]= 0.;
	          scal_prod_(param_int[nd], ipt_ind_dm_thread, krylov_krp1, krylov_i, ipt_norm_kry[ithread-1]);
                  printf("normvect   %g %d  \n", ipt_norm_kry[ithread-1] , i);
               }*/
               nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
      }
lhs_end = omp_get_wtime();
if (ithread==1) printf("Normalisation_vect     %g \n", lhs_end-lhs_beg);
lhs_beg = omp_get_wtime();

/*
    // norm L2 pour verrif bug
    for (E_Int th = 0; th < Nbre_thread_actif; th++) { normL2_sum += ipt_norm_kry[th];}
    normL2_sum = sqrt(normL2_sum);
    printf("normvect   %f %d  \n",normL2_sum , ithread);
*/

    kr++;

//#include  "Compute/Linear_solver/verif_vectkrylov.cpp"

    continue_gmres = (K_FUNC::E_abs(ipt_VectG[kr]) / save) > epsi_linear;
    //#include "Compute/Linear_solver/verif_vectkrylov.cpp"

  }//loop kr
   // fin loop vecteur krylov
   //  
   //
   //


  //Resolution de Y par remontee
  E_Float* Hessenberg   = ipt_Hessenberg;
  for (E_Int i = kr - 1; i >= 0; i--)
    {
       E_Int iline   = i *(num_max_vect - 1);
       value = 0.;
       ipt_norm_kry[ithread-1]= 0.;
       //#pragma omp simd for
       #pragma omp for
       for (E_Int j = kr - 1; j > i; j--) {  ipt_norm_kry[ithread-1] -= ipt_VectY[j] * ipt_Hessenberg[ iline + j ]; }

       #pragma omp single
       {
        /* for (E_Int i = 0; i < kr ; i++) */
        /*   { */
        /* 	E_Float* Hessenberg_i   = ipt_Hessenberg + i * (num_max_vect - 1); */
        /* 	printf("hess apr \n " ); */
        /* 	for (E_Int j = 0; j < kr - 1; j++) { printf(" %f %d %d ", Hessenberg_i[j], i,j ); } */
        /* 	printf("\n" ); */
        /*   } */


        for (E_Int th = 0; th < Nbre_thread_actif; th++) { value += ipt_norm_kry[th];}
        ipt_VectY[i] = (ipt_VectG[i] + value) / Hessenberg[ iline + i ];

        //printf("VecY  %g %d  \n", ipt_VectY[i], i );
        } //end single

   }//end loop i
lhs_end = omp_get_wtime();
if (ithread==1) printf("Remontee               %g \n", lhs_end-lhs_beg);
lhs_beg = omp_get_wtime();

  if (ithread ==1) printf("Residu GMRES =  %14.12g, target= %14.12g,  Nit_kry= %d, Nb_restart= %d \n",K_FUNC::E_abs(ipt_VectG[ kr ]) / save, epsi_linear, kr, restart );

  //continue_gmres = K_FUNC::E_abs(ipt_VectG[kr]) > epsi_linear;
  //cout << "continue_gmres = " << continue_gmres << endl;

  kr++;
  //Calcul de X solution du GMRES, stockage dans drodm
  shift_zone =0; nd_current =0;
  for (E_Int nd = 0; nd < nidom; nd++)
    {
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
       // 
       // inverser les loop et vectoriser avec simd reduction
       // 
       prod_mat_vect_(param_int[nd], ipt_ind_dm_thread, iptkrylov[nd],
                      ipt_VectY, iptdrodm + shift_zone, kr);
       nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
     shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
    }//loop zone
lhs_end = omp_get_wtime();
if (ithread==1) printf("prod_mat_vect          %g \n", lhs_end-lhs_beg);
lhs_beg = omp_get_wtime();

      shift_coe =0; shift_zone =0; nd_current =0;
      for (E_Int nd  = 0; nd < nidom; nd++)
        {
#include "HPC_LAYER/OMP_MODE_BEGIN.h"
#include "Compute/LU/prep_lussor.h"
	  indice = nd * nb_subzone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + ithread - 1;

            E_Float* krylov_in = iptdrodm + shift_zone;
	    E_Float* krylov_out= iptdrodm + shift_zone;
	    mjrnewton = 1;

	    invlu_(nd                     , nitcfg      ,nitrun, param_int[nd], param_real[nd],
	    	   ipt_ind_dm_thread      , ipt_ind_dm_thread       , mjrnewton             ,
	    	   iptrotmp[nd]           , iptro_ssiter[nd]        , krylov_in             , krylov_out            ,
	    	   ipti[nd]               , iptj[nd]                , iptk[nd]              ,
	    	   iptventi[nd]           , iptventj[nd]            , iptventk[nd]          ,
	    	   iptcoe  + shift_coe    , ipt_ssor_shift          , iptssortmp[nd]        , ipt_ssor_size[ indice ]);

            nd_current +=1;
#include "HPC_LAYER/OMP_MODE_END.h"
          shift_zone = shift_zone + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ];
          shift_coe  = shift_coe  + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ_COE ];
         }//loop zone
lhs_end = omp_get_wtime();
if (ithread==1) printf("invlu final            %g \n", lhs_end-lhs_beg);
lhs_beg = omp_get_wtime();
  //
  //
  //init_gramm schmidt again normL2
  //
  //
  //
  //Include de print de verif

//#include "verif_Ax-b.cpp"
