         E_Int sizefluD= recv_fluxRc[irac][1];
         E_Int nobcR   = recv_fluxRc[irac][3];

         E_Int pt_bcs   = param_int[ NoR ][PT_BC];
         E_Int nb_bc    = param_int[ NoR ][pt_bcs];
         E_Int adrFlu   = param_int[ NoR ][pt_bcs + 1 + nobcR + nb_bc];
         E_Float* fluxR =param_real[ NoR ] + adrFlu;
         
         E_Int sizefluR= sizefluD/4;
         if(param_int[ NoR ][NIJK+4]==0) { sizefluR= sizefluD/2;} // Cas 2D

         E_Int shift   = nvars_loc * (nbRcvPts + shift_fluR);
         shift_fluR   += sizefluR;
         //printf("No raccord Recep %d , nbRcvPts: %d ,  sizefluR: %d, shift_fluR: %d \n", nobcR, nbRcvPts, sizefluR, shift_fluR);
         #pragma omp for nowait
         for (E_Int l = 0; l < sizefluR*nvars_loc; l++)
           {
            fluxR[ l ] = recv_frp[irac][ l + shift];
            //if(NoR==0) {printf("flux recep %f,  %d %d \n", fluxR[ l ], l, l/sizefluR );}
           }

