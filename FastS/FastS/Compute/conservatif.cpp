   E_Int nitcfg_loc = 1;

    //E_Int nidom_cons = 147;
    //E_Int nidom_cons = 136;
    //E_Int nidom_cons = 132;
    //E_Int nidom_cons = 273;
    E_Int nidom_cons = 274;
    //E_Int nidom_cons =10;
    //for (E_Int nd = 0; nd < nidom_cons; nd++)

    debit[ithread-1                      ]=0;
    debit[ithread-1 + Nbre_thread_actif  ]=0;
    debit[ithread-1 + Nbre_thread_actif*2]=0;
    debit[ithread-1 + Nbre_thread_actif*3]=0;
    debit[ithread-1 + Nbre_thread_actif*4]=0;
    debit[ithread-1 + Nbre_thread_actif*5]=0;
    debit[ithread-1 + Nbre_thread_actif*6]=0;
    debit[ithread-1 + Nbre_thread_actif*7]=0;
    debit[ithread-1 + Nbre_thread_actif*8]=0;
    debit[ithread-1 + Nbre_thread_actif*9]=0;
    debit[ithread-1 + Nbre_thread_actif*10]=0;
    debit[ithread-1 + Nbre_thread_actif*11]=0;

    nd_current=0;

    for (E_Int nd = 0; nd < nidom; nd++)
            {

             if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
             {
               E_Float* ipt_CL = iptro_CL[nd];

               E_Int shift_deb = Nbre_thread_actif*6;
               if(nd >= nidom_cons) shift_deb= 0;
               E_Float* ipt_debit = debit+ ithread-1 + shift_deb;

                      E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;   //nidom_loc(nssiter)
                      E_Int  nb_subzone    = ipt_nidom_loc [nitcfg_loc-1];                                           //nbre sous-zone a la sousiter courante
                      for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
                        {
                         E_Int ndo   = nd;

                         E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg_loc-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
 
                         E_Int* ipt_ind_dm_thread;
                         if (omp_mode == 1)
                         { 
                           E_Int       Ptomp = param_int[nd][PT_OMP];
                           E_Int  PtrIterOmp = param_int[nd][Ptomp +nitcfg_loc -1];   
                           E_Int  PtZoneomp  = param_int[nd][PtrIterOmp + nd_subzone];

                           Nbre_thread_actif_loc = param_int[nd][ PtZoneomp  + Nbre_thread_actif ];
                           ithread_loc           = param_int[nd][ PtZoneomp  +  ithread -1       ] +1 ;
                           ipt_ind_dm_thread     = param_int[nd] + PtZoneomp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;
  
                           if (ithread_loc == -1) { nd_current++; continue;}
                         }
                         else
                         { 
                           E_Int* ipt_topology_socket = ipt_topology       + (ithread-1)*3;
                           E_Int* ipt_ind_dm_socket   = ipt_ind_dm_omp     + (ithread-1)*12;
                           ipt_ind_dm_thread   = ipt_ind_dm_socket  +6;

                           E_Int lmin = 10;
                           if (param_int[nd][ITYPCP] == 2) lmin = 4;

                           indice_boucle_lu_(ndo, ithread, Nbre_thread_actif, lmin,
                                             ipt_ind_dm_loc,
                                             ipt_topology_socket, ipt_ind_dm_thread);
                         }// omp_mode

               //if(ithread==1)printf("cococo %d \n",nd);
               cp_conservatif_(param_int[nd], param_real[nd], Nbre_thread_actif, ipt_ind_dm_thread, ipt_CL, iptvol[nd], iptCellN[nd] ,ipt_debit);

               nd_current +=1;
                     
#include       "HPC_LAYER/OMP_MODE_END.h"

             }//maillage structure

            }//fin boucle zone


#pragma omp barrier
#pragma omp single
   {
       for (E_Int ne= 0; ne <  6 ; ne++){ 
        for (E_Int i = 1;  i < Nbre_thread_actif ; i++) 
        { debit[ 0                   +ne*Nbre_thread_actif ]  += debit[i                       +ne*Nbre_thread_actif ];
          debit[ Nbre_thread_actif*6 +ne*Nbre_thread_actif ]  += debit[i + Nbre_thread_actif*6 +ne*Nbre_thread_actif ]; }}

    masse[0]=  debit[ Nbre_thread_actif*6 ];
    masse[1]= param_real[0][PSIROE+1];
    masse[2]= param_real[0][PSIROE+2];
    varMasse =  0.5*(3.*masse[0]-4*masse[1]+masse[2])/ param_real[0][DTC];
    ro_corr =  0.666666666666666*param_real[0][DTC]/debit[Nbre_thread_actif *11]*(-varMasse +  param_real[0][EPSI_INFLOW]);

    if(nitrun < 3){ro_corr=0;}

    for (E_Int nd1= 0;          nd1 <  nidom_cons   ; nd1++){  param_real[nd1][ROTATION] =  ro_corr;}
    for (E_Int nd1= nidom_cons; nd1 <  nidom ; nd1++){         param_real[nd1][ROTATION] =  0;}

    if(nitcfg==30)
     {  param_real[0][PSIROE+2] =masse[1];
        param_real[0][PSIROE+1] =masse[0];
     }

   }//omp single

    if (ithread==1 && nitcfg==30) printf("ro_corr   = %d  %.15f  %.15f  %.15f  %.15f \n",nitcfg, ro_corr,-varMasse +  param_real[0][EPSI_INFLOW] ,varMasse,  param_real[0][EPSI_INFLOW] );

    /*
    nd_current=0;
    for (E_Int nd = 0; nd < nidom_cons; nd++)
      {

       if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
        {
         E_Float* ipt_CL = iptro_CL[nd];

         E_Int shift_deb = Nbre_thread_actif*6;
         if(nd >= nidom_cons) shift_deb= 0;
         E_Float* ipt_debit = debit+ ithread-1 + shift_deb;

#include   "HPC_LAYER/OMP_MODE_NITLOC_BEGIN.h"

            if(nitcfg!=31) {corr_conservatif_(param_int[nd], param_real[nd], Nbre_thread_actif, ipt_ind_dm_thread, ipt_CL, iptvol[nd], iptCellN[nd] ,ipt_debit);}

            nd_current +=1;
                     
#include    "HPC_LAYER/OMP_MODE_END.h"
         }//maillage structure
       }//fin boucle zone
    */


//Verif correction
/*
    debit[ithread-1                      ]=0;
    debit[ithread-1 + Nbre_thread_actif  ]=0;
    debit[ithread-1 + Nbre_thread_actif*2]=0;
    debit[ithread-1 + Nbre_thread_actif*3]=0;
    debit[ithread-1 + Nbre_thread_actif*4]=0;
    debit[ithread-1 + Nbre_thread_actif*5]=0;
    debit[ithread-1 + Nbre_thread_actif*6]=0;
    debit[ithread-1 + Nbre_thread_actif*7]=0;
    debit[ithread-1 + Nbre_thread_actif*8]=0;
    debit[ithread-1 + Nbre_thread_actif*9]=0;
    debit[ithread-1 + Nbre_thread_actif*10]=0;
    debit[ithread-1 + Nbre_thread_actif*11]=0;

    nd_current=0;

    for (E_Int nd = 0; nd < nidom; nd++)
            {

             if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
             {
               E_Float* ipt_CL = iptro_CL[nd];

               E_Int shift_deb = Nbre_thread_actif*6;
               if(nd >= 251) shift_deb= 0;
               E_Float* ipt_debit = debit+ ithread-1 + shift_deb;

#include       "HPC_LAYER/OMP_MODE_NITLOC_BEGIN.h"

               //if(ithread==1)printf("cococo \n");
               cp_conservatif_(param_int[nd], param_real[nd], Nbre_thread_actif, ipt_ind_dm_thread, ipt_CL, iptvol[nd], iptCellN[nd] ,ipt_debit);

               nd_current +=1;
                     
#include       "HPC_LAYER/OMP_MODE_END.h"

             }//maillage structure

            }//fin boucle zone


#pragma omp barrier
#pragma omp single
   {
       for (E_Int ne= 0; ne <  6 ; ne++){ 
        for (E_Int i = 1;  i < Nbre_thread_actif ; i++) 
        { debit[ 0                   +ne*Nbre_thread_actif ]  += debit[i                       +ne*Nbre_thread_actif ];
          debit[ Nbre_thread_actif*6 +ne*Nbre_thread_actif ]  += debit[i + Nbre_thread_actif*6 +ne*Nbre_thread_actif ]; }}

    masse[0]=  debit[ Nbre_thread_actif*6 ];
    masse[1]= param_real[0][PSIROE+1];
    masse[2]= param_real[0][PSIROE+2];
    varMasse =  0.5*(3.*masse[0]-4*masse[1]+masse[2])/ param_real[0][DTC];
    ro_corr =  0.6666*param_real[0][DTC]/debit[Nbre_thread_actif *11]*(-varMasse +  param_real[0][EPSI_INFLOW]);
        

    if(nitcfg==30)
     {  param_real[0][PSIROE+2] =masse[1];
        param_real[0][PSIROE+1] =masse[0];
     }
   }//omp single

    if (ithread==1) printf("vr_corr    = %d  %.15f  %.15f  %.15f  %.15f \n",nitcfg, ro_corr,-varMasse +  param_real[0][EPSI_INFLOW] ,varMasse,  param_real[0][EPSI_INFLOW] );
    //if (ithread==1) printf("vr_corr    = %d  %.15f  %.15f  %.15f  %.15f %.15f \n",nitcfg, masse[0],masse[1],masse[2] ,varMasse,  param_real[0][EPSI_INFLOW] );
 */
