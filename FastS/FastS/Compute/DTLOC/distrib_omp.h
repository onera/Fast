      if(omp_mode==1)
      {
        E_Int       Ptomp = param_int[NoD][PT_OMP];
        E_Int  PtrIterOmp = param_int[NoD][Ptomp +nstep -1];   
        E_Int  PtZoneomp  = param_int[NoD][PtrIterOmp + 0];

        E_Int count[nb_socket];
        for (E_Int c = 0; c < nb_socket; c++){count[c]=0;}
        //on determine sur quelle socket est la zone 
        for (E_Int th = 0; th < __NUMTHREADS__; th++)
           { 
            E_Int socket = th/core_per_socket;
            E_Int th_loc = param_int[NoD][ PtZoneomp  +  th      ] +1 ;
            if(th_loc!=-1) count[socket]+=1;
           } 

        // si 2 socket, la zone est calcule en paralle par tous les coeur du socket
        if  (nb_socket==2) 
             { 
               E_Int th_min=1; E_Int th_max = activ_core_per_socket[0];  E_Int socket_tg = 0;
               if (count[1] > count[0]) {th_min = activ_core_per_socket[0] +1; th_max =  activ_core_per_socket[0] +  activ_core_per_socket[1]; socket_tg=1;}

               if ( th_min <= ithread && th_max >= ithread ) {Nbre_thread_actif_loc = activ_core_per_socket[socket_tg]; ithread_loc = ithread- socket_tg*activ_core_per_socket[0];}
               else{  continue;}
               //if(ithread==1 && irac==0)  printf("verif soc_tg= %d, ith_loc= %d %d , Nthread_actif= %d , NoZone= %d %d %d %d \n",  socket_tg, ithread_loc, ithread, Nbre_thread_actif_loc, NoD, count[0], count[1], nstep);
             }

        else { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread; }
      }
