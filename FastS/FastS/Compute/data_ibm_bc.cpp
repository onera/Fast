    //flux conservatif: calcul moyenne vitesse interface
    E_Int pt_flu = param_int[nd][IBC_PT_FLUX];
    E_Int shift =6*Nfamily +1 ;
    E_Int ithread =1;
    E_Int Nthread_max=1;

    if(nd ==0)
     {
        for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
        {
          E_Int shift_fam = fam*7*Nthread_max;  E_Int shift_th = (ithread-1)*7;  E_Int shift = shift_fam + shift_th;
          for ( E_Int i  = 0; i <  7; i++ ) {flux [shift +i]=0; }
        }
     }

    /*
    if (pt_flu != -1)
     {
      for ( E_Int fam  = 0; fam <  Nfamily; fam++ )
       {
        E_Float* iptflu = flux + (ithread-1)*7 + fam*Nthread_max*7;
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
               cp_debit_ibm_(nd, idir, neq_mtr, ithread, Nthread_max, nitcfg, param_int[nd],  param_real[nd],
                             size_fen, facelist, iptro_CL[nd], iptijk, iptvol[nd], iptflu);

              shift += size_fen;
             }
         }//dir
       }//family
      }
      */

     ///mise a jour moyenne plan Lund si necessaire
     E_Int pt_bcs = param_int[nd][PT_BC];
     E_Int nb_bc  = param_int[nd][ pt_bcs ];
     E_Float* ipt_data;
     for ( E_Int ndf = 0; ndf < nb_bc; ndf++ )
     {
        E_Int pt_bc  = param_int[nd][pt_bcs+ 1 + ndf];

        E_Int idir   = param_int[nd][pt_bc + BC_IDIR];
        E_Int nbdata = param_int[nd][pt_bc + BC_NBDATA];
        E_Int bc_type= param_int[nd][pt_bc + BC_TYPE];

        E_Int* iptsize_data = param_int[nd] + pt_bc + BC_NBDATA + 1;
        E_Int* ind_fen      = param_int[nd] + pt_bc + BC_FEN;
        E_Int  inc_bc[3];
        if ( idir <= 2 ) 
        { inc_bc[0] = ind_fen[3] - ind_fen[2] + 1; // nombre element de la fenetre dans la direction J
          inc_bc[1] = ind_fen[2]; // debut indice j
          inc_bc[2] = ind_fen[4]; // debut indice k
        }  
        else if ( idir <= 4 )
        { inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
          ind_fen[0];             // debut indice i
          inc_bc[2] = ind_fen[4]; // debut indice k
        }  
        else
        {inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
         inc_bc[1] = ind_fen[0]; // debut indice i
         inc_bc[2] = ind_fen[2]; // debut indice j
        }

        if ( nbdata != 0 ) ipt_data = param_real[nd] + param_int[nd][pt_bcs + 1 + ndf + nb_bc];

        E_Int* ipt_ind_dm_loc  = ipt_ind_dm[nd]  + (nitcfg-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;
        //Lund
        if(bc_type==19 and nitcfg==nitcfg_last-1) 
         { 
            E_Float* iptAvgPlanLund = ipt_data + 5*iptsize_data[0];
            E_Float* iptParamLund   = ipt_data + 10*iptsize_data[0];

            mj_lund_planrecyl_(nd, idir, param_int[nd] , ind_fen, ipt_ind_dm_loc, inc_bc, iptsize_data[0],
                            iptParamLund,  iptro[nd], iptAvgPlanLund); 
         }
        //wallmodel
        else if(bc_type==31)
         {  E_Float* iptAvgPlan; E_Float* iptParam_wmles;
            E_Int sz_data;
            if(iptsize_data[0]==1){iptAvgPlan     = ipt_data + iptsize_data[0]; 
                                   iptParam_wmles = ipt_data + iptsize_data[1];
                                   sz_data =(iptsize_data[1]-iptsize_data[0])/5/iptParam_wmles[0];}
            else{iptAvgPlan = ipt_data;  iptParam_wmles = ipt_data + iptsize_data[0]; sz_data= iptsize_data[0]/5/iptParam_wmles[0]; }

            //No du snap a sauvegarder
            iptParam_wmles[2]+=1;
            if(iptParam_wmles[2] == iptParam_wmles[0]) {iptParam_wmles[2]=1;}

            E_Float* iptSnapPlan = iptAvgPlan + 5*sz_data*int(iptParam_wmles[2]);
            //E_Float* iptSnapPlan = iptAvgPlan + 5*sz_data*int(2); // modif vmoy

            // codage sequentiel!!!!!!
            //if(nd==1) printf("No snap/ nitrun %f %f  %d nitcfg= %d \n", iptParam_wmles[2], iptParam_wmles[2], nitrun, nitcfg ) ;
            E_Int lprint =0;
            if(nitcfg==nitcfg_last-1 && nitrun%50==0) lprint =1;

            mj_wallmodel_plan_(nd, idir, param_int[nd] , ind_fen, ipt_ind_dm_loc, inc_bc, sz_data, lprint,
                               iptParam_wmles,  iptro_CL[nd], iptAvgPlan, iptSnapPlan );
         }//wallmodel
     }//ndf
