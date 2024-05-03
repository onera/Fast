    if(nitcfg<=30)
     {
        //printf("Masse    = %d  %.12f   %.12f  %.12f \n",nitrun, debit[0   ]+debit[ 24*6 ],debit[0   ] ,debit[ 24*6 ]);
        //printf("MomentumX= %d  %.12f   %.12f  %.12f \n",nitrun, debit[24  ]+debit[ 24*7 ],debit[24  ] ,debit[ 24*7 ]);
        //printf("MomentumY= %d  %.12f   %.12f  %.12f \n",nitrun, debit[24*2]+debit[ 24*8 ],debit[24*2] ,debit[ 24*8 ]);
        //printf("Energy   = %d  %.12f   %.12f  %.12f \n",nitrun, debit[24*4]+debit[ 24*10],debit[24*4] ,debit[ 24*10]);

       E_Float   c4 =  5. / 6.;
       E_Float   c5 =  2. / 6.;
       E_Float   c6 = -1. / 6.;
       E_Float vparoi_old[Nfamily*2];
       E_Float vparoi[Nfamily*2];
       E_Float vparoi_tot=0; E_Float vparoi_tot1=0; E_Float vparoi_tot2=0;
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
       //for ( E_Int fam  = 0; fam <  1; fam++ )
         {
          //E_Float* iptflu = flux + (ithread-1)*10 + fam*Nthread_max*7;
          E_Float* iptflu = flux + fam*neqFlu;
          vparoi_old[fam        ]= (c5*iptflu[0] +c4*iptflu[1] + c6*iptflu[2]);
          vparoi_old[fam+Nfamily]= (c4*iptflu[0] +c5*iptflu[1] + c6*iptflu[6]);
          vparoi_tot  +=vparoi_old[fam];
          vparoi_tot1 +=vparoi_old[fam+Nfamily];
          //printf("flux paroiD0  %.15f %.15f fam= %d , nstep/nitrun: %d %d \n", vparoi_old[fam], vparoi_tot, fam, nitcfg, nitrun);
          //printf("flux paroiG0  %.15f  %.15f  \n", vparoi_old[fam+Nfamily], vparoi_tot1);
         }
        //printf("flux paroiD0  %.15f %.15f %.15f %.15f %.15f \n", vparoi_old[0], vparoi[1],vparoi[2],vparoi[3], vparoi_tot);
        //printf("flux paroiG0  %.15f %.15f %.15f %.15f %.15f \n", vparoi_old[0+Nfamily], vparoi[1+Nfamily],vparoi[2+Nfamily],vparoi[3+Nfamily], vparoi_tot1);
        //printf("flux paroiD0  %.15f %.15f %.15f %d %d  \n", vparoi[0], vparoi[1], vparoi_tot, nitcfg, nitrun);
        //printf("flux paroiG0  %.15f %.15f %.15f  \n", vparoi[0+Nfamily], vparoi[1+Nfamily], vparoi_tot1);
       /*
       */ 

       //estimation perte paroi apres 1er correction de musker
       E_Int ithread =1;
       E_Int Nthread_max=1;
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
        {
         E_Int shift_fam = fam*neqFlu*Nthread_max;  E_Int shift_th = (ithread-1)*neqFlu;  E_Int shift = shift_fam + shift_th;
         for ( E_Int i  = 0; i < neqFlu ; i++ ) {flux [shift +i   ]=0; }
        }

       for (E_Int nd = 0;  nd < nidom ; nd++) 
       {
       E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
       E_Int shift =6*Nfamily +1 ;
       if (pt_flu != -1)
         {
         vparoi_tot=0;  vparoi_tot1=0; vparoi_tot2=0;
         for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
         //for ( E_Int fam  = 0; fam <  1; fam++ )
          {
           E_Float* iptflu = flux + (ithread-1)*neqFlu + fam*Nthread_max*neqFlu;
           E_Float* iptamor= amor + fam*2;
           for ( E_Int idir = 1; idir<= 6; idir++ )
            {
                 // printf("verif %d , fam= %d , idir= %d , size= %d %d \n", pt_flu,fam ,idir,param_int[nd][pt_flu+idir+fam*6],nitcfg);
                 //printf("verif flux777 %.15f %d %d %d %d \n", iptflu1[6], nd, fam, idir, (ithread-1)*neqFlu + 3*Nthread_max*neqFlu);
          
              E_Int size_fen =  param_int[nd][pt_flu+idir +fam*6];    
              if (size_fen != 0 )
                {
                  //printf("verif %d , fam= %d , idir= %d , size= %d %d \n", pt_flu,fam ,idir,param_int[nd][pt_flu+idir+fam*6],nitcfg);
                  E_Int* facelist = param_int[nd] + pt_flu  + shift;
                  E_Float* iptijk; E_Float* ipventijk; E_Int neq_mtr;
                  if      ( idir <= 2 ) { iptijk    = ipti[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventi[nd]; } 
                  else if ( idir <= 4 ) { iptijk    = iptj[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventj[nd]; }
                  else                  { iptijk    = iptk[nd]; neq_mtr   = param_int[nd][NEQ_K ]; ipventijk = iptventk[nd]; }
                  cp_debit_ibm_(nd, idir, neq_mtr, ithread, Nthread_max, nitcfg, param_int[nd],  param_real[nd],
                                size_fen, facelist, iptro_CL[nd], iptijk, iptvol[nd], iptflu);

                 shift += size_fen;
                }
            }//dir

          //if(nitcfg<=30 and nd == nidom-1) {printf("flux family   %.15f  %.15f %.15f %.15f %.15f %.15f %.15f %.15f %d %d  \n", iptflu[5], iptflu[10], iptflu[11],  amor[fam*2  ], amor[fam*2 +1 ], iptflu[7], iptflu[8], iptflu[3], fam, nitcfg);}

           vparoi[fam]        = (c5*iptflu[0]  +c4*iptflu[1] + c6*iptflu[2]);
           vparoi[fam+Nfamily]= (c4*iptflu[0]  +c5*iptflu[1] + c6*iptflu[6]);
           vparoi_tot        += iptflu[7];
           vparoi_tot1       += iptflu[8];
           vparoi_tot2       += iptflu[10];

           iptamor[0] =   (-iptflu[11]*0.5-iptflu[8])/iptflu[12]/c5;                     // 4*residu/c5/(rhol+rhor)
           iptamor[1] = ( (-iptflu[11]*0.5-iptflu[7])/iptflu[12] -c4*iptamor[0])/c6;

           //amor[fam*2  ] = vparoi[fam]/(c5*iptflu[3]);
           //amor[fam*2+1] = (vparoi[fam+Nfamily] -c4/c5*vparoi[fam])/(c6*iptflu[3]);

          }//family
         }//flu
       //if(nitcfg==30  and nd == nidom-1) {printf("flux familes: rho Ul/r   %.15f %.15f , flux num %.15f %.15f %d  \n", vparoi_tot, vparoi_tot1, vparoi_tot+vparoi_tot1 ,vparoi_tot2,  nitrun);}
       }//nidom

//for (E_Int loop = 0; loop < 3; loop++)
for (E_Int loop = 0; loop < 0; loop++)
     {

#pragma omp parallel default(shared)
 {
#ifdef _OPENMP
   E_Int  ithread           = omp_get_thread_num() +1;
   E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
   E_Int ithread = 1;
   E_Int Nbre_thread_actif = 1;
#endif
   #pragma omp for
   for (E_Int nd1 = 0; nd1 < nidom_ibc; nd1++)
     { 
       E_Int nd =  nd_ibc[nd1];
       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
       E_Int  nb_subzone    = ipt_nidom_loc [nitcfg_stk-1];                                       //nbre sous-zone a la sousiter courante

       if (autorisation_bc[nd] == 1 && nb_subzone > 0)
        {
           E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
           E_Int shift =6*Nfamily +1 ;
           if (pt_flu != -1)
            {
             E_Int ithread_local = 1; E_Int Nbre_thread_actif_local = 1;
             for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
         //for ( E_Int fam  = 0; fam <  1; fam++ )
              {
                E_Float* iptamor= amor + fam*2;
                E_Float* iptflu = flux + fam*neqFlu;
                for ( E_Int idir = 1; idir<= 6; idir++ )
                {
                  E_Int size_fen =  param_int[nd][pt_flu+idir +fam*6];    
                  if (size_fen != 0 )
                   {
                     E_Int* facelist = param_int[nd] + pt_flu  + shift;
                     E_Float* iptijk; E_Float* ipventijk; E_Int neq_mtr;
                     if      ( idir <= 2 ) { iptijk    = ipti[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventi[nd]; } 
                     else if ( idir <= 4 ) { iptijk    = iptj[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventj[nd]; }
                     else                  { iptijk    = iptk[nd]; neq_mtr   = param_int[nd][NEQ_K ]; ipventijk = iptventk[nd]; }

                     E_Int ipass= 0;
                     corr_debit_ibm_(nd, idir, neq_mtr, ithread_local, Nbre_thread_actif_local,  param_int[nd], param_real[nd], ipass,
                                    fam, nitcfg, iptamor, size_fen, facelist, iptro_CL[nd], iptijk, iptvol[nd], iptCellN[nd], iptflu);

                     shift += size_fen;
                   }
                }//dir
               }//fam
             }
        }//autorisation
     }//loop zone
 }//fin zone omp

       ithread =1;
       Nthread_max=1;
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
        {
         E_Int shift_fam = fam*neqFlu*Nthread_max;  E_Int shift_th = (ithread-1)*neqFlu;  E_Int shift = shift_fam + shift_th;
         for ( E_Int i  = 0; i <  neqFlu; i++ ) {flux [shift +i   ]=0; }
        }


       for (E_Int nd = 0;  nd < nidom ; nd++) 
       {
       E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
       E_Int shift =6*Nfamily +1 ;
       if (pt_flu != -1)
         {
         for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
         //for ( E_Int fam  = 0; fam <  1; fam++ )
          {
           E_Float* iptflu = flux + (ithread-1)*neqFlu + fam*Nthread_max*neqFlu;
           for ( E_Int idir = 1; idir<= 6; idir++ )
            {
                 // printf("verif %d , fam= %d , idir= %d , size= %d %d \n", pt_flu,fam ,idir,param_int[nd][pt_flu+idir+fam*6],nitcfg);
                 //printf("verif flux777 %.15f %d %d %d %d \n", iptflu1[6], nd, fam, idir, (ithread-1)*neqFlu + 3*Nthread_max*neqFlu);
          
              E_Int size_fen =  param_int[nd][pt_flu+idir +fam*6];    
              if (size_fen != 0 )
                {
                  //printf("verif %d , fam= %d , idir= %d , size= %d %d \n", pt_flu,fam ,idir,param_int[nd][pt_flu+idir+fam*6],nitcfg);
                  E_Int* facelist = param_int[nd] + pt_flu  + shift;
                  E_Float* iptijk; E_Float* ipventijk; E_Int neq_mtr;
                  if      ( idir <= 2 ) { iptijk    = ipti[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventi[nd]; } 
                  else if ( idir <= 4 ) { iptijk    = iptj[nd]; neq_mtr   = param_int[nd][NEQ_IJ]; ipventijk = iptventj[nd]; }
                  else                  { iptijk    = iptk[nd]; neq_mtr   = param_int[nd][NEQ_K ]; ipventijk = iptventk[nd]; }
                  cp_debit_ibm_(nd, idir, neq_mtr, ithread, Nthread_max, nitcfg, param_int[nd],  param_real[nd],
                                size_fen, facelist, iptro_CL[nd], iptijk, iptvol[nd], iptflu);

                 shift += size_fen;

          //printf("debit: nd= %d  %.12f %.12f %.12f %.12f %d \n",nd, flux[0],flux[1],flux[2],flux[3], facelist[0]  );
                }
            }//dir
          //if(nitcfg<=30 and nd == nidom-1) {printf("flux fam loo  %.15f  %.15f %.15f %.15f %.15f %.15f %.15f %.15f %d %d  \n", iptflu[5], iptflu[10], iptflu[11],  amor[fam*2  ], amor[fam*2 +1 ], iptflu[7], iptflu[8], iptflu[3], fam, loop);}
          }//family
         }//flu
       }//nidom



       //E_Float   c4 =  5. / 6.;
       //E_Float   c5 =  2. / 6.;
       //E_Float   c6 = -1. / 6.;

       vparoi_tot=0;  vparoi_tot1=0;
       for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
       //for ( E_Int fam  = 0; fam <  1; fam++ )
         {
          //E_Float* iptflu = flux + (ithread-1)*10 + fam*Nthread_max*neqFlu;
          E_Float* iptflu = flux + fam*neqFlu;
          E_Float* iptamor= amor + fam*2;

          //vparoi[fam]        = (c5*iptflu[0]  +c4*iptflu[1] + c6*iptflu[2]);
          //vparoi[fam+Nfamily]= (c4*iptflu[0]  +c5*iptflu[1] + c6*iptflu[6]);
          vparoi[fam]  = iptflu[10]; 
          vparoi_tot   +=vparoi[fam];
          vparoi_tot1  +=iptflu[13];

          //if(nitcfg==30) {printf("flux paroiD1  %.15f %.15f fam= %d , amor= %.15f %.15f , loop/nstep/nitrun: %d %d %d \n", vparoi[fam],vparoi[fam+Nfamily] , fam, amor[fam*2],amor[fam*2+1],loop, nitcfg, nitrun); }
          if(nitcfg==30 and loop == 2 and fam == 7) {printf("flux family loop %d  %.15f  %.15f \n",fam,  vparoi_tot, vparoi_tot1);}
          //printf("flux verif  %.15f  %.15f %.15f \n", vparoi_old[fam        ], -c4*iptflu[1] - c6*iptflu[2],(c5*iptflu[0]) );

          iptamor[0] =   (-iptflu[11]*0.5-iptflu[8])/iptflu[12]/c5;                     // 4*residu/c5/(rhol+rhor)
          iptamor[1] = ( (-iptflu[11]*0.5-iptflu[7])/iptflu[12] -c4*iptamor[0])/c6;

          E_Float amorT =    (-iptflu[17]*0.5-iptflu[16])/iptflu[14]/c5;
          E_Float amorT1 =(  (-iptflu[17]*0.5-iptflu[15])/iptflu[14] -c4*amorT)/c6;

          //if(nitcfg<=30 and loop == 2 ) {printf("flux family loop %d  %.15f  %.15f %.15f  %.15f \n",fam, iptamor[0], iptamor[1], amorT, amorT1 );}
          //if(nitcfg<=30 and loop == 2 ) {printf("flux family loop %d  %.15f  %.15f %.15f   \n",fam, iptamor[0], iptamor[1], (-iptflu[17]*0.5-iptflu[16])/iptflu[14]/c5 );}
          //if(nitcfg<=30 and loop == 2 ) {printf("flux family loop %d  %.15f  %.15f %.15f  %.15f \n",fam, iptamor[0], iptamor[1], amorT, amorT1 );}
          //amor[fam*2  ] = vparoi[fam]/(c5*iptflu[3]);
          //amor[fam*2+1] = (vparoi[fam+Nfamily] -c4/c5*vparoi[fam])/(c6*iptflu[3]);

         }
        //printf("flux paroiD1  %.15f %.15f  %d %d  \n",  vparoi_tot, vparoi_tot1, nitcfg, nitrun);
        //printf("flux paroiD1  %.15f %.15f %.15f %d %d  \n", vparoi[0], vparoi[1], vparoi_tot, nitcfg, nitrun);
        //printf("flux paroiG1  %.15f %.15f %.15f  \n", vparoi[0+Nfamily], vparoi[1+Nfamily], vparoi_tot1);
        //printf("flux paroiD1  %.15f %.15f  %.15f %d %d \n", vparoi[0], vparoi[1], vparoi_tot, nitcfg, lexit_lu);
        //printf("flux paroiG1  %.15f %.15f  %.15f \n", vparoi[0+Nfamily], vparoi[1+Nfamily], vparoi_tot1);

  } //loop

     }//nit=30
