         //verrou rhs
         E_Int type = 2;
         //si calcul residu on verifie que toutes les souszones sont pretes
         if(   nb_subzone > 1 && (param_int[nd][ ITYPCP] != 2 || param_int[nd][ DTLOC ]== 1) )
          {
            E_Int size = param_int[nd][NEQ]*param_int[nd][NDIMDX]; flush_real_( size , iptdrodm + shift_zone);

            for (E_Int ntask_loc = 0; ntask_loc < nbtask; ntask_loc++)
             {
              E_Int pttask_loc   = ptiter + ntask_loc*(6+Nbre_thread_actif*7);

              E_Int nd_loc                 = ipt_omp[ pttask_loc                         ];
              E_Int Nbre_thread_actif_loc1 = ipt_omp[ pttask_loc + 2 + Nbre_thread_actif ];

              if(nd == nd_loc)
               { for (E_Int th = 0; th < Nbre_thread_actif_loc1; th++) 
                 { E_Int* verrou_lhs_thread= verrou_lhs + ntask_loc*Nbre_thread_actif + th; verrou_c_( verrou_lhs_thread, type); } 
               }
             }
          }
         else  //on verifie uniquememnt la souszone courante: ntask
          { if(ithread_loc != -1)
             { for (E_Int th = 0; th < Nbre_thread_actif_loc; th++) 
                 { E_Int* verrou_lhs_thread= verrou_lhs + ntask*Nbre_thread_actif + th; verrou_c_( verrou_lhs_thread, type); }
             }
          }
