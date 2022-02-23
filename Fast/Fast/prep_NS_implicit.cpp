//
//mise a jour Nombre sous_iter pour implicit (simplification gestion OMP)
//
for (E_Int nd = 0; nd < nidom; nd++)
  {
    if( (param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1 || lssiter_verif ==1 ) && param_int[nd][ ITYPZONE] != 4)
     {
      E_Int* ipt_nisdom_residu   =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter;                //nisdom_residu(nssiter)
      E_Int* ipt_it_bloc         =  ipt_ind_dm[nd]      + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter*2;    //it_bloc(nidom)

      if(ipt_nisdom_residu[nitcfg-1] != 0) ipt_it_bloc[0] +=1;
     }
   }

//
//Calcul taille tableau ssor par thread 
//
E_Int nbtask = ipt_omp[nitcfg-1]; 
E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

for (E_Int ntask = 0; ntask < nbtask; ntask++)
{
#ifdef _OPENMP
E_Int  Nbre_thread_actif = __NUMTHREADS__;
#else
E_Int Nbre_thread_actif = 1;
#endif
  E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
  E_Int nd = ipt_omp[ pttask ];

  if( (param_int[ nd ][ NB_RELAX ] > 1 || param_int[ nd ][ LU_MATCH ]==1) && param_int[nd][ ITYPZONE] != 4)
   {
    E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)

    E_Int* ipt_inddm_omp;
    E_Int nd_subzone      = ipt_omp[ pttask + 1 ];
    E_Int Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
    for (E_Int i = 0; i < Nbre_thread_actif; i++)
      {
       E_Int ithread_loc     = ipt_omp[ pttask + 2 + i] +1 ;
       ipt_inddm_omp         = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

       E_Int indice = nd * mx_sszone * Nbre_thread_actif + nd_subzone * Nbre_thread_actif + i;
       if (ithread_loc != -1)
         { ipt_ssor_size[indice]=(ipt_inddm_omp[1] - ipt_inddm_omp[0] +1 +2*param_int[nd][ NIJK +3 ])*
                                 (ipt_inddm_omp[3] - ipt_inddm_omp[2] +1 +2*param_int[nd][ NIJK +3 ])*
                                 (ipt_inddm_omp[5] - ipt_inddm_omp[4] +1 +2*param_int[nd][ NIJK +4 ]);
         }
       else { ipt_ssor_size[indice]=0; }
      } // loop threads
   } // if relax
} // loop task
