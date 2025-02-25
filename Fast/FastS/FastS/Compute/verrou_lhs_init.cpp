   for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
       E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
       ithread_loc  = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;

       if (ithread_loc == -1) {continue;}

       //Reinitialisation verrou omp rhs
       E_Int l =  ntask*mx_synchro*Nbre_thread_actif  + (ithread_loc-1)*mx_synchro;
       for (E_Int i = 0;  i < mx_synchro ; i++) { ipt_lok[ l + i ]  = 0; }

       E_Int type = 4;
       //Init verrou rhs pour chaque sous zone et chaque thread actif:  init val to zero
       E_Int* verrou_lhs_thread = verrou_lhs +  ntask*Nbre_thread_actif + ithread_loc -1; 
       verrou_c_( verrou_lhs_thread, type);

       //Reinitialisation verrou omp pour calcul residu avant lhs
       verrou_lhs_thread        = verrou_lhs + (nbtask + ntask)*Nbre_thread_actif + ithread_loc -1; 
       verrou_c_( verrou_lhs_thread, type );
     }
