    if(nitcfg==30)
     {
        printf("Masse    = %d  %.12f   %.12f  %.12f \n",nitrun, debit[0   ]+debit[ 24*6 ],debit[0   ] ,debit[ 24*6 ]);
        //printf("MomentumX= %d  %.12f   %.12f  %.12f \n",nitrun, debit[24  ]+debit[ 24*7 ],debit[24  ] ,debit[ 24*7 ]);
        //printf("MomentumY= %d  %.12f   %.12f  %.12f \n",nitrun, debit[24*2]+debit[ 24*8 ],debit[24*2] ,debit[ 24*8 ]);
        //printf("Energy   = %d  %.12f   %.12f  %.12f \n",nitrun, debit[24*4]+debit[ 24*10],debit[24*4] ,debit[ 24*10]);

       E_Float   c4 =  5. / 6.;
       E_Float   c5 =  2. / 6.;
       E_Float   c6 = -1. / 6.;
       E_Float vparoi[Nfamily*2];
       E_Float vparoi_tot=0; E_Float vparoi_tot1=0;
       /*
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
         {
          //E_Float* iptflu = flux + (ithread-1)*10 + fam*Nthread_max*7;
          E_Float* iptflu = flux + fam*7;
          vparoi[fam        ]= (c5*iptflu[0] +c4*iptflu[1] + c6*iptflu[2]);
          vparoi[fam+Nfamily]= (c4*iptflu[0] +c5*iptflu[1] + c6*iptflu[6]);
          vparoi_tot  +=vparoi[fam];
          vparoi_tot1 +=vparoi[fam+Nfamily];
         }
        printf("flux paroiD0  %.15f %.15f %.15f %.15f %.15f \n", vparoi[0], vparoi[1],vparoi[2],vparoi[3], vparoi_tot);
        printf("flux paroiG0  %.15f %.15f %.15f %.15f %.15f \n", vparoi[0+Nfamily], vparoi[1+Nfamily],vparoi[2+Nfamily],vparoi[3+Nfamily], vparoi_tot1);
       */ 
       E_Int ithread =1;
       E_Int Nthread_max=1;
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
        {
         E_Int shift_fam = fam*7*Nthread_max;  E_Int shift_th = (ithread-1)*7;  E_Int shift = shift_fam + shift_th;
         for ( E_Int i  = 0; i <  7; i++ ) {flux [shift +i   ]=0; }
        }


       for (E_Int nd = 0;  nd < nidom ; nd++) 
       {
       E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
       E_Int shift =6*Nfamily +1 ;
       if (pt_flu != -1)
         {
         for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
          {
           E_Float* iptflu = flux + (ithread-1)*7 + fam*Nthread_max*7;
           for ( E_Int idir = 1; idir<= 6; idir++ )
            {
                 // printf("verif %d %d %d %d %d \n", pt_flu,fam ,idir,param_int[nd][pt_flu+idir+fam*6],nitcfg);
                  //printf("verif flux777 %.15f %d %d %d %d \n", iptflu1[6], nd, fam, idir, (ithread-1)*7 + 3*Nthread_max*7);
          
              E_Int size_fen =  param_int[nd][pt_flu+idir +fam*6];    
              if (size_fen != 0 )
                {
                  E_Int* facelist = param_int[nd] + pt_flu  + shift;
                  E_Float* iptijk; E_Float* ipventijk; E_Int neq_mtr;
                  if      ( idir <= 2 ) { iptijk    = ipti[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventi[nd]; } 
                  else if ( idir <= 4 ) { iptijk    = iptj[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventj[nd]; }
                  else                  { iptijk    = iptk[nd]; neq_mtr   = param_int[nd][NEQ_K ]; ipventijk = iptventk[nd]; }
                  cp_debit_ibm_(nd, idir, neq_mtr, ithread, Nthread_max, nitcfg, param_int[nd],  param_real[nd],
                                size_fen, facelist, iptro_CL[nd], iptijk, iptvol[nd], iptflu);

                 shift += size_fen;

          //if(nd == 42 )printf("debit: nit/nd= %d  %d  %d %d %.12f %.12f %d \n",idir, Nfamily, fam, size_fen, flux[4],flux[5], facelist[0]  );
                }
            }//dir
          }//family
         }//flu
       }//nidom

       vparoi_tot=0;  vparoi_tot1=0;
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
         {
          //E_Float* iptflu = flux + (ithread-1)*10 + fam*Nthread_max*7;
          E_Float* iptflu = flux + fam*7;
          vparoi[fam]= (c5*iptflu[0] +c4*iptflu[1] + c6*iptflu[2]);
          vparoi[fam+Nfamily]= (c4*iptflu[0] +c5*iptflu[1] + c6*iptflu[6]);
          vparoi_tot  +=vparoi[fam];
          vparoi_tot1 +=vparoi[fam+Nfamily];
         }
        printf("flux paroiD1  %.15f %.15f %.15f %.15f %.15f \n", vparoi[0], vparoi[1],vparoi[2],vparoi[3], vparoi_tot);
        printf("flux paroiG1  %.15f %.15f %.15f %.15f %.15f \n", vparoi[0+Nfamily], vparoi[1+Nfamily],vparoi[2+Nfamily],vparoi[3+Nfamily], vparoi_tot1);

     }//nit=30
