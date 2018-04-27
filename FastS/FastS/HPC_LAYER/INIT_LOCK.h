   //E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
   E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if( Nbre_thread_actif < Nbre_socket ) Nbre_socket = 1;

   E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
   E_Int socket          = (ithread-1)/thread_parsock +1;
 //E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

        for (E_Int nd = 0; nd < nidom; nd++)
          {
             E_Int l =  nd*mx_synchro*Nbre_thread_actif  + (ithread-1)*mx_synchro;
             for (E_Int i = 0;  i < mx_synchro ; i++)
                { ipt_lok[ l + i ]  = 0; }
          }
          #pragma omp barrier
