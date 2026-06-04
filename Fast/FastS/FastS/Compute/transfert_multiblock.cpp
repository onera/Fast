E_Int cycl;
E_Float deb_dtlocal; E_Float tmps;

E_Int flag_passage2= 0;
E_Int numpassage   = 1;
E_Int rk           = param_int[0][RK];
E_Int exploc       = param_int[0][EXPLOC];

  if (param_int[0][EXPLOC] == 0) //   dt constant
  {
      for (E_Int nopass = 0; nopass < iptdtloc[12]; nopass++)
        {
         //printf("N0pass trans %d \n",  nopass+1);fflush(0);
         K_FASTC::setInterpTransfersFast(iptro_CL, vartype, int_tc[nopass], real_tc[nopass] , param_int, param_real, ipt_omp,
                                     linelets_int, linelets_real, it_target, nidom, ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage, nopass );
        }
  }
  else if (param_int[0][EXPLOC] == 1)     // explicit local instationnaire
  {
     for (E_Int nopass = 0; nopass < iptdtloc[12]; nopass++)
       {
         K_FASTS::dtlocal2para_c(iptro, iptrotmp, int_tc[nopass], param_int, param_real, iptdrodm, iptcoe, stock, drodmstock, constk, nitcfg, taille_tabs, nidom);
       }

     for (E_Int nd=0; nd < nidom; nd++)
       {
        cycl = nssiter/param_int[nd][LEVEL];
        if (nitcfg%cycl == 0 and nitcfg != nssiter)
	   {
	     E_Float* ptsave  = iptro[nd]; 
	     iptro[nd] = iptrotmp[nd]; 
	     iptrotmp[nd] = ptsave;
	   } 

	 if (nitcfg%cycl == cycl/2 and (nitcfg/cycl)%2==1 ) { flag_passage2=1; }
       }

     if (nitcfg != nssiter)
      {  		
        for (E_Int nopass = 0; nopass < iptdtloc[12]; nopass++)
         {
           K_FASTC::setInterpTransfersFast(iptro_CL, vartype, int_tc[nopass], real_tc[nopass] , param_int, param_real, ipt_omp,
                                        linelets_int, linelets_real, it_target, nidom , ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage, nopass);
         }
      }

     for (E_Int nopass = 0; nopass < iptdtloc[12]; nopass++)
       {
         K_FASTS::recup3para_c(iptro, int_tc[nopass], param_int, stock, nitcfg, taille_tabs, nidom);
       }

     //BC_local(iptro, iptrotmp, param_int_tc, param_real_tc, param_int, param_real, iptdrodm, iptcoe, ipt_ind_CL, ipti, iptj, iptk, iptx, ipty, iptz, iptventi, iptventj, iptventk, stock, drodmstock, constk, nitcfg, omp_mode, taille_tabs, nidom);

     numpassage=2;
     if (flag_passage2==1)
      {
        for (E_Int nopass = 0; nopass < iptdtloc[12]; nopass++)
         {
           K_FASTC::setInterpTransfersFast(iptro_CL, vartype, int_tc[nopass], real_tc[nopass] , param_int, param_real, ipt_omp,
                                        linelets_int, linelets_real, it_target, nidom , ipt_timecount, mpi, nitcfg, nssiter, rk, exploc, numpassage, nopass);
         }
        flag_passage2=0;
      }   
   } // fin boucle test dtlocal

  else{ printf("Erreur transfert_multibloc.cpp: exploc soucis \n"); exit(0);}
