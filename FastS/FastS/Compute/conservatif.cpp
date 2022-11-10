   E_Int nitcfg_loc = 1;

    //E_Int nidom_cons = 147;
    //E_Int nidom_cons = 136;
    //E_Int nidom_cons = 132;
    //E_Int nidom_cons = 273;
    //E_Int nidom_cons = 274;
    E_Int nidom_cons = nidom;
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

    nbtask = ipt_omp[nitcfg_loc-1]; 
    ptiter = ipt_omp[nssiter+ nitcfg_loc-1];

    for (E_Int ntask = 0; ntask < nbtask; ntask++)
      {
         E_Int pttask     = ptiter + ntask*(6+Nbre_thread_actif*7);
         E_Int nd         = ipt_omp[ pttask ];
         E_Int nd_subzone = ipt_omp[ pttask + 1 ];

         shift_zone=0; shift_coe=0;
         for (E_Int n = 0; n < nd; n++)
          {
           shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
           shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
          }

         if (param_int[nd][ITYPZONE] != 4 and param_int[nd][IFLOW] != 4)  //on skippe les eventuelles zone non structurees ou LBM
           {
             E_Float* ipt_CL = iptro_CL[nd];

             E_Int shift_deb = Nbre_thread_actif*6;
             if(nd >= nidom_cons) shift_deb= 0;
             E_Float* ipt_debit = debit+ ithread-1 + shift_deb;


             E_Int* ipt_topo_omp; E_Int* ipt_ind_dm_thread;
 
             ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
             Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
             ipt_ind_dm_thread     = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

             if (ithread_loc == -1) {continue;}

               //if(ithread==1)printf("cococo %d \n",nd);
               cp_conservatif_(param_int[nd], param_real[nd], Nbre_thread_actif, ipt_ind_dm_thread, ipt_CL, iptvol[nd], iptCellN[nd] ,ipt_debit);

           }//maillage structure

      }//fin boucle task  

   E_Int  it_print =30;
   if(param_int[0][ITYPCP]==2){it_print =3;}
#pragma omp barrier
#pragma omp single
   {
       for (E_Int ne= 0; ne <  6 ; ne++){ 
        for (E_Int i = 1;  i < Nbre_thread_actif ; i++) 
        { debit[ 0                   +ne*Nbre_thread_actif ]  += debit[i                       +ne*Nbre_thread_actif ];
          debit[ Nbre_thread_actif*6 +ne*Nbre_thread_actif ]  += debit[i + Nbre_thread_actif*6 +ne*Nbre_thread_actif ]; }}

          //printf("veriffff %f %f \n", debit[ 0], debit[ Nbre_thread_actif*6]);

    masse[0]=  debit[ Nbre_thread_actif*6 ];
    masse[1]= param_real[0][PSIROE+1];
    masse[2]= param_real[0][PSIROE+2];
    varMasse =  0.5*(3.*masse[0]-4*masse[1]+masse[2])/ param_real[0][DTC];
    //ro_corr =  0.666666666666666*param_real[0][DTC]/debit[Nbre_thread_actif *11]*(-varMasse +  param_real[0][EPSI_INFLOW]);
    ro_corr =  0.666666666666666*param_real[0][DTC]/debit[Nbre_thread_actif *11]*(-varMasse);

    if(nitrun < 3){ro_corr=0;}

    //for (E_Int nd1= 0;          nd1 <  nidom_cons   ; nd1++){  param_real[nd1][ROTATION] =  ro_corr;}
    //for (E_Int nd1= nidom_cons; nd1 <  nidom ; nd1++){         param_real[nd1][ROTATION] =  0;}

    if(nitcfg==it_print)
     {  param_real[0][PSIROE+2] =masse[1];
        param_real[0][PSIROE+1] =masse[0];
     }

   }//omp single

    //if (ithread==1 && nitcfg==30) printf("ro_corr   = %d  %.15f  %.15f  %.15f  %.15f \n",nitcfg, ro_corr,-varMasse +  param_real[0][EPSI_INFLOW] ,varMasse,  param_real[0][EPSI_INFLOW] );
    if (ithread==1 && nitcfg==it_print) printf("ro_corr   = %d  %.15f  %.15f  %.15f  %.15f \n",nitrun, ro_corr,-varMasse +  param_real[0][EPSI_INFLOW] ,varMasse/masse[0],masse[0]  );

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
