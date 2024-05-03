   //E_Int Nbre_socket   = NBR_SOCKET;                       
   E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
   if (Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

   E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
   E_Int socket          = (ithread-1)/thread_parsock +1;
   //E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;
