//barrier pour attendre toutes les zones
if (Nfamily != 0)
{ 
#pragma omp barrier
}
     //estimation perte paroi apres 1er correction de musker
     //E_Int ithread =1;
     //E_Int Nthread_max=1;
     for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
        {
         E_Int shift_fam = fam*neqFlu*Nthread_max;  E_Int shift_th = (ithread-1)*neqFlu;  E_Int shift = shift_fam + shift_th;
         for ( E_Int i  = 0; i < neqFlu ; i++ ) {flux [shift +i   ]=0; }
        }

     E_Int nitcfg_loc =1;
     nbtask = ipt_omp[nitcfg_loc-1]; 
     ptiter = ipt_omp[nssiter+ nitcfg_loc-1];

     for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
       E_Int pttask          = ptiter + ntask*(6+Nbre_thread_actif*7);
       E_Int nd              = ipt_omp[ pttask ];
       ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
       Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
       if (ithread_loc == -1) {continue;}

       shift_zone=0; shift_coe=0;
       for (E_Int n = 0; n < nd; n++)
         {
            shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
            shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
         }

       E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
       E_Int shift =6*Nfamily +1 ;
       if (pt_flu != -1)
         {
         for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
          {
           E_Float* iptflu = flux + (ithread-1)*neqFlu + fam*Nthread_max*neqFlu;

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
                  //cp_debit_ibm_(nd, idir, neq_mtr, ithread_loc, Nbre_thread_actif_loc , nitcfg, param_int[nd],  param_real[nd],
                  //              size_fen, facelist, iptro_CL[nd], iptijk, iptvol[nd], iptflu);
                  cp_corr_debit_ibm_(nd, idir, neq_mtr, ithread_loc, Nbre_thread_actif_loc , nitcfg, param_int[nd],  param_real[nd],
                                size_fen, facelist, iptro_CL[nd], iptijk, iptcoe  + shift_coe , iptdrodm  + shift_zone,  iptflu);

                 shift += size_fen;
                }
            }//dir

          }//family
         }//flu
       }//ntask  

if (Nfamily != 0)
{ 
#pragma omp barrier
}
/*
// barrier pour attendre toutes les zones
#pragma omp barrier
#pragma omp single
   {
     for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
        {
         E_Int shift_fam = fam*neqFlu*Nthread_max;  
         for ( E_Int i  = 0; i < neqFlu ; i++ )
           {
            for ( E_Int ith  = 1; ith < Nthread_max ; ith++ )
              {
               E_Int shift_th  = ith*neqFlu; 
               flux [shift_fam +i   ] +=flux [shift_fam + shift_th +i ]; 
              }//Nthreads
           }//neq
        } //fam
   }//single

  E_Float tx=0; 
  E_Float ty=0; 
  E_Float tz=0; 
  for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
  //for ( E_Int fam  = 0; fam <  1; fam++ )
        {
          E_Int q =  fam*neqFlu*Nthread_max + 19;
          tx += flux[q];
          ty += flux[q+1];
          tz += flux[q+2];
          if(ithread==1){ printf("flux %d %d %.15f %.15f %.15f \n",  fam, nitcfg, flux[q], flux[q+1],  flux[q+2] ); }
        }
     if(ithread==1){ printf("bilan normal %d %.15f %.15f %.15f \n",  nitcfg, tx, ty, tz); }
*/

/*
     //correction bilan 
     for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
       E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
       E_Int nd     = ipt_omp[ pttask ];
       ithread_loc  = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
       Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
       if (ithread_loc == -1) {continue;}

       shift_zone=0; shift_coe=0;
       for (E_Int n = 0; n < nd; n++)
         {
            shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
            shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
         }

       E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
       E_Int shift =6*Nfamily +1 ;
       if (pt_flu != -1)
         {
         for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
         //for ( E_Int fam  = 0; fam <  1; fam++ )
          {
           //E_Float* iptflu = flux + (ithread-1)*neqFlu + fam*Nthread_max*neqFlu;
           E_Float* iptflu = flux +  fam*Nthread_max*neqFlu;
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
                  corr_bilan_ibm_(nd, idir, neq_mtr, ithread_loc, Nbre_thread_actif_loc, nitcfg, param_int[nd],  param_real[nd],
                                size_fen, facelist, iptdrodm + shift_zone , iptijk, iptcoe  + shift_coe, iptflu);

                 shift += size_fen;
                }
            }//dir

          }//family
         }//flu
       }//ntask 
*/
