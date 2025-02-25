         pos = param_int_tc[shift_rac + nrac*17];
         E_Int* fluD = param_int_tc + pos + 4*nbflu;
         E_Int idir    = fluD[0];
         E_Int sizefluD= fluD[1];
         E_Int nobcD   = fluD[2];
         E_Int nobcR   = fluD[3];

         E_Int pt_bcs= param_int[ NoD ][PT_BC];
         E_Int nb_bc = param_int[ NoD ][pt_bcs];
         E_Int pt_bc = param_int[ NoD ][pt_bcs + 1 + nobcD];
         E_Int adrFlu= param_int[ NoD ][pt_bcs + 1 + nobcD + nb_bc];
         E_Int* fen  = param_int[ NoD ] + pt_bc + BC_FEN;
         E_Float* fluxD = param_real[ NoD ] + adrFlu;

         pt_bcs= param_int[ NoR ][PT_BC];
         nb_bc = param_int[ NoR ][pt_bcs];
         adrFlu= param_int[ NoR ][pt_bcs + 1 + nobcR + nb_bc];
         E_Float* fluxR = param_real[ NoR ] + adrFlu;
         //E_Int sizefluR= sizefluD/4;
         // Cas 2D
         if(param_int[ NoR ][NIJK+4]==0)
           { 
            E_Int sizefluR= sizefluD/2;
            if(idir<=2)
              {
               E_Int jmax= (fen[3]-fen[2]+1)/2;
               #pragma omp for nowait
               for (E_Int j = 0; j < jmax; j++)
                {
                  E_Int lr  = j;
                  E_Int ld0 = j*2;
                  E_Int ld1 = j*2 +1;
                 //if(NoR==48) {printf("jr %d, kr %d, lr %d , ld0 %d , ld1 %d \n", j, k, lr, ld0, ld1);}
                  fluxR[lr            ] = fluxD[ld0            ];
                  fluxR[lr            ]+= fluxD[ld1            ];
                  fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                  fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                  fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                  fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                  fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                  fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                  fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                  fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
                  fluxR[lr +sizefluR*5] = fluxD[ld0 +sizefluD*5];
                  fluxR[lr +sizefluR*5]+= fluxD[ld1 +sizefluD*5];
                  //printf("fluR %f %f %f %f %f \n", fluxR[lr],fluxR[lr +sizefluR],fluxR[lr +sizefluR*2], fluxR[lr +sizefluR*3],fluxR[lr +sizefluR*4]);
                }
              }
            else if(idir<=4)
              {
               E_Int imax= (fen[1]-fen[0]+1)/2;
               #pragma omp for nowait
               for (E_Int i = 0; i < imax; i++)
                {
                  E_Int lr  = i;
                  E_Int ld0 = i*2;
                  E_Int ld1 = i*2 +1;
                
                  fluxR[lr            ] = fluxD[ld0            ];
                  fluxR[lr            ]+= fluxD[ld1            ];
                  fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                  fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                  fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                  fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                  fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                  fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                  fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                  fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
                  fluxR[lr +sizefluR*5] = fluxD[ld0 +sizefluD*5];
                  fluxR[lr +sizefluR*5]+= fluxD[ld1 +sizefluD*5];
                }
              }
           }
         else  //3D
           { 
            E_Int sizefluR= sizefluD/4;
            if(idir<=2)
             {
              E_Int jmax= (fen[3]-fen[2]+1)/2;
              E_Int kmax= (fen[5]-fen[4]+1)/2;
              E_Int inck = jmax*2;
              #pragma omp for collapse(2) nowait
              for (E_Int k = 0; k < kmax; k++)
               {
               for (E_Int j = 0; j < jmax; j++)
                {
                 E_Int lr  = j   + k*jmax;
                 E_Int ld0 = j*2 + k*2*inck;
                 E_Int ld1 = ld0 +1 ;
                 E_Int ld2 = ld0 +inck;
                 E_Int ld3 = ld1 +inck;

                 fluxR[lr            ] = fluxD[ld0];
                 fluxR[lr            ]+= fluxD[ld1];
                 fluxR[lr            ]+= fluxD[ld2];
                 fluxR[lr            ]+= fluxD[ld3];
                 fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld2 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld3 +sizefluD  ];
                 fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld2 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld3 +sizefluD*2];
                 fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld2 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld3 +sizefluD*3];
                 fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld2 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld3 +sizefluD*4];
                 fluxR[lr +sizefluR*5] = fluxD[ld0 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld1 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld2 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld3 +sizefluD*5];
                }
               }
             }
            else if(idir<=4)
             {
              E_Int imax= (fen[1]-fen[0]+1)/2;
              E_Int kmax= (fen[5]-fen[4]+1)/2;
              E_Int inck = imax*2;
              #pragma omp for collapse(2) nowait
              for (E_Int k = 0; k < kmax; k++)
               {
               for (E_Int i = 0; i < imax; i++)
                {
                 E_Int lr  = i   + k*imax;
                 E_Int ld0 = i*2 + k*2*inck;
                 E_Int ld1 = ld0 +1 ;
                 E_Int ld2 = ld0 +inck;
                 E_Int ld3 = ld1 +inck;
                 fluxR[lr            ] = fluxD[ld0];
                 fluxR[lr            ]+= fluxD[ld1];
                 fluxR[lr            ]+= fluxD[ld2];
                 fluxR[lr            ]+= fluxD[ld3];
                 fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld2 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld3 +sizefluD  ];
                 fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld2 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld3 +sizefluD*2];
                 fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld2 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld3 +sizefluD*3];
                 fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld2 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld3 +sizefluD*4];
                 fluxR[lr +sizefluR*5] = fluxD[ld0 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld1 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld2 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld3 +sizefluD*5];
                }
               }
             }
            else 
             {
              E_Int imax= (fen[1]-fen[0]+1)/2;
              E_Int jmax= (fen[3]-fen[2]+1)/2;
              E_Int incj = imax*2;
              #pragma omp for collapse(2) nowait
              for (E_Int j = 0; j < jmax; j++)
               {
               for (E_Int i = 0; i < imax; i++)
                {
                 E_Int lr  = i   + j*imax;
                 E_Int ld0 = i*2 + j*2*incj;
                 E_Int ld1 = ld0 +1 ;
                 E_Int ld2 = ld0 +incj;
                 E_Int ld3 = ld1 +incj;

                 fluxR[lr            ] = fluxD[ld0];
                 fluxR[lr            ]+= fluxD[ld1];
                 fluxR[lr            ]+= fluxD[ld2];
                 fluxR[lr            ]+= fluxD[ld3];
                 fluxR[lr +sizefluR  ] = fluxD[ld0 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld1 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld2 +sizefluD  ];
                 fluxR[lr +sizefluR  ]+= fluxD[ld3 +sizefluD  ];
                 fluxR[lr +sizefluR*2] = fluxD[ld0 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld1 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld2 +sizefluD*2];
                 fluxR[lr +sizefluR*2]+= fluxD[ld3 +sizefluD*2];
                 fluxR[lr +sizefluR*3] = fluxD[ld0 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld1 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld2 +sizefluD*3];
                 fluxR[lr +sizefluR*3]+= fluxD[ld3 +sizefluD*3];
                 fluxR[lr +sizefluR*4] = fluxD[ld0 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld1 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld2 +sizefluD*4];
                 fluxR[lr +sizefluR*4]+= fluxD[ld3 +sizefluD*4];
                 fluxR[lr +sizefluR*5] = fluxD[ld0 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld1 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld2 +sizefluD*5];
                 fluxR[lr +sizefluR*5]+= fluxD[ld3 +sizefluD*5];
                }
               }
             }
           }//2d/3d

