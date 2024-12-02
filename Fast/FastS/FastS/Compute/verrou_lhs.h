          if(lssiter_verif ==1   && ( param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
           {
             E_Int type = 2;
             for (E_Int ntask_loc = 0; ntask_loc < nbtask; ntask_loc++)
              {
                E_Int pttask_loc     = ptiter + ntask_loc*(6+Nbre_thread_actif*7);
                E_Int nd_loc         = ipt_omp[ pttask_loc ];
                E_Int nd_subzone_loc = ipt_omp[ pttask_loc + 1 ];
          
                E_Int Nbre_thread_actif_loc1 = ipt_omp[ pttask_loc + 2 + Nbre_thread_actif ];
                E_Int ithread_loc1           = ipt_omp[ pttask_loc + 2 + ithread -1 ] +1 ;

                if(nd == nd_loc &&  nd_subzone_loc == 0 && barrier_residu==0 && (nb_subzone > 1 || ithread_loc1 != -1 ) )
                {
                    if (ithread_loc1 == -1) {ithread_loc1 =1;} // astuce pour verrou pour arreter thread souszone multiple. verrou unique pour souzone=0

                    for (E_Int th = 0; th < Nbre_thread_actif_loc1; th++) 
                     { 
                       E_Int shift = (ithread_loc1-1 + th)%Nbre_thread_actif_loc1;
       
                       E_Int* verrou_lhs_thread= verrou_lhs + (nbtask + ntask_loc)*Nbre_thread_actif + shift; 
                       verrou_c_( verrou_lhs_thread, type); 
                     }
                }
              }//loop taskloc
           } //sinon residu pas bon en omp_mode=1

          //Wait thread avant attaquer LU 
          if(lssiter_verif ==0 && kimpli==1)
           {
              E_Int type = 2;
              for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
                { E_Int* verrou_lhs_thread= verrou_lhs + ntask*Nbre_thread_actif + th; verrou_c_( verrou_lhs_thread, type); }
           }
