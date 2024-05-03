#ifdef NB_SOCKET
  E_Int nb_socket=NB_SOCKET;
#else
  E_Int nb_socket=1;
#endif
#ifdef CORE_PER_SOCK
  E_Int core_per_socket=CORE_PER_SOCK;
#else
  E_Int core_per_socket=__NUMTHREADS__;
#endif

//On dertermine le nombre de socket reelememnt utilisees
nb_socket=__NUMTHREADS__/core_per_socket;
if(nb_socket==0) nb_socket +=1;
 
E_Float cout_moyen = 0; 
for (E_Int so = 0; so < nb_socket; so++) //loop sur nb socket
  { 
    E_Float cout_moyen = 0;
    E_Int th_persock =core_per_socket;
    E_Int th_deb = so*core_per_socket; 
    E_Int th_fin = (so+1)*core_per_socket; 
    if(th_fin > __NUMTHREADS__) th_fin = __NUMTHREADS__;

    for (E_Int th = th_deb; th < th_fin; th++) // loop sur thread du socket pour calculer cout moyen sur les ssiter
      { 
       //if(nitcfg==nssiter-1)
       if(nitcfg==nitcfg_last)
         { E_Float cout1 = 0; E_Float cout2 = 0;
           E_Int nit_moy = nitcfg-1;
           if (param_int[0][ITYPCP]==2) nit_moy = nssiter;
           for (E_Int it = 1; it <= nit_moy; it++)
             { E_Int ptTimer = th*2 + (it-1)*threadmax_sdm*2;
                      cout1 += timer_omp[ptTimer   ];
                      cout2 += timer_omp[ptTimer +1];
             }
             //printf("Thread=  %d, CPU rhs= %f, CPU lhs= %f, CPU tot= %f, nitcfg= %d, nitrun= %d \n", th+1, cout1,cout2, cout1 +cout2, nitcfg, nitrun); 
         }
        E_Int ptTimer = th*2 + (nitcfg-1)*threadmax_sdm;
        cout_moyen       += timer_omp[ptTimer ]+timer_omp[ptTimer +1];
      }


    E_Float marge =1.010;
    cout_moyen   /= float(th_fin -th_deb +1);

    E_Int Nbre_thread_actif = threadmax_sdm;

    for (E_Int ithread = th_deb+1; ithread <= th_fin; ithread++)
    {
     E_Int ptTimer = (ithread-1)*2 + (nitcfg-1)*Nbre_thread_actif*2;
     E_Float cout_th = timer_omp[ptTimer ]+timer_omp[ptTimer +1];
     //si thread lent, on cherche a lui enlever du travail
     if (cout_th > marge*cout_moyen) 
     { 
        //Recherche thread rapide de la socket
        E_Float cout_min = 100000000.;
        E_Int th_min;
        for (E_Int th = th_deb+1 ; th <= th_fin; th++)
          {  
              E_Int ptTimer = (th-1)*2 + (nitcfg-1)*Nbre_thread_actif*2;
              cout_th       = timer_omp[ptTimer]+timer_omp[ptTimer+1];
              if (cout_th < cout_min){ th_min = th;  cout_min = cout_th;}
          } 

        //estimation cout total ithread
        E_Int cells_tg=0; E_Int count_zone=0;
        for (E_Int nd = 0; nd < nidom; nd++)
        { 
         E_Int * ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
         E_Int nb_subzone    = ipt_nidom_loc [nitcfg-1];      

         for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
         {
            E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp;
            if (omp_mode == 1)
            { 
              E_Int       Ptomp = param_int[nd][PT_OMP];
              E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
              E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

              E_Int ithread_loc = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
              ipt_inddm_omp     = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;

              if (ithread_loc != -1)
                  {E_Int c1    = (ipt_inddm_omp[1]-ipt_inddm_omp[0]+1)*(ipt_inddm_omp[3]-ipt_inddm_omp[2]+1)*(ipt_inddm_omp[5]-ipt_inddm_omp[4]+1);
                   cells_tg   += c1;
                   count_zone +=1;
                   //printf("th %d : size zone= %d , nijk= %d %d %d  nd= %d, cell_tg= %d \n", ithread, c1, ipt_inddm_omp[1],ipt_inddm_omp[3],ipt_inddm_omp[5], nd, cells_tg);
                   //if(ithread==24) printf("th 24: size zone= %d , nijk= %d %d %d  nd= %d, cell_tg= %d \n", c1, ipt_inddm_omp[1],ipt_inddm_omp[3],ipt_inddm_omp[5], nd, cells_tg);
                  }
            }
         }
        } 
        E_Float cout_per_cell = cout_th /float(cells_tg);

        //estimation zone a enlever; on force le choix dans une zone mono threader pour le moment
        E_Int nd_tg;  E_Int nd_subzone_tg; E_Float dist =10000000.; E_Int cell_ajoutee; E_Int search =0;
        for (E_Int nd = 0; nd < nidom; nd++)
        { 
         E_Int * ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
         E_Int nb_subzone    = ipt_nidom_loc [nitcfg-1];      
         for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
         {
            E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp;
            if (omp_mode == 1)
            { 
              E_Int       Ptomp = param_int[nd][PT_OMP];
              E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg -1];   
              E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

              E_Int Nbre_thread_actif_loc = param_int[nd][ PtZoneomp  + Nbre_thread_actif ];

              E_Int ithread_loc     = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
              ipt_inddm_omp         = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;
             
              if (ithread_loc != -1 && Nbre_thread_actif_loc ==1) 
              {  
                 E_Int cells = (ipt_inddm_omp[1]-ipt_inddm_omp[0]+1)*(ipt_inddm_omp[3]-ipt_inddm_omp[2]+1)*(ipt_inddm_omp[5]-ipt_inddm_omp[4]+1);
                 E_Float cell_cible = (cout_th -cout_moyen)*cout_per_cell;

                 E_Float dist_loc = ( cells- cell_cible)*( cells- cell_cible);
                 if(dist_loc < dist) {  dist = dist_loc; nd_subzone_tg= nd_subzone; nd_tg=nd; cell_ajoutee=cells; search=1;}
              }
            }
         }
        } 

        if (omp_mode == 1 && search ==1)
        { 
              E_Int       Ptomp = param_int[nd_tg][PT_OMP];
              E_Int  PtrIterOmp = param_int[nd_tg][Ptomp +nitcfg -1];   
              E_Int  PtZoneomp  = param_int[nd_tg][PtrIterOmp + nd_subzone_tg];

              E_Int ithread_loc                            = param_int[nd_tg][ PtZoneomp  +  ithread -1       ] +1 ;
              param_int[nd_tg][ PtZoneomp  +  ithread -1 ] = -2;
              param_int[nd_tg][ PtZoneomp  +  th_min  -1 ] = ithread_loc-1; // a blinder si thread deja actif sur cette zone
               
              E_Int ptTimer = (th_min-1)*2 + (nitcfg-1)*Nbre_thread_actif*2;

              timer_omp[ptTimer ] += cell_ajoutee*cout_per_cell;
      
          //printf("thread lent = %d, thread rapide= %d, zone= %d, subs= %d, Nb_subs= %d, Ptr %d \n",  ithread, th_min, nd_tg,nd_subzone_tg, count_zone, PtZoneomp  +  ithread -1);
        }

   }//  test cout_th
 }//loop thread


        }//loop socket
        if(nitcfg==nitcfg_last)
        {

         for (E_Int it = 1; it <= nssiter; it++)
          {
          for (E_Int th = 0; th < threadmax_sdm; th++) 
           { E_Int ptTimer = th*2 + (it-1)*threadmax_sdm*2; 
             timer_omp[ptTimer   ]=0;
             timer_omp[ptTimer+1 ]=0;
             //printf("mise zero = %d  %d  %d  %d \n", nssiter, it, nitcfg, ptTimer );
           }//loop thread
          }//loop it
        }

