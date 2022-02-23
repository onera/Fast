        #pragma omp barrier
        for (E_Int nd = 0; nd < nidom; nd++)
          {
             E_Int l =  nd*mx_synchro*Nbre_thread_actif  + (ithread-1)*mx_synchro;
             for (E_Int i = 0;  i < mx_synchro ; i++)
                { ipt_lok[ l + i ]  = 0; }
          }
