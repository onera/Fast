E_Float meax       = 0;
E_Float meay       = 0;
E_Float meaz       = 0;

#pragma omp parallel default(shared) reduction(+:meax,meay,meaz)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
    E_Float rhs_begin        = omp_get_wtime();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
    E_Float rhs_begin       = 0;
#endif
    //E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
    E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
    if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

    E_Int Nbre_thread_actif_loc, ithread_loc;
    if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
    else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

    E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
    E_Int socket          = (ithread-1)/thread_parsock +1;
    E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

    E_Int* ipt_topology_socket    = ipt_topology       + (ithread-1)*3;
    E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm       + (ithread-1)*3;
    E_Int* ipt_ind_CL_thread      = ipt_ind_CL         + (ithread-1)*6;
    E_Int* ipt_ind_CL119          = ipt_ind_CL         + (ithread-1)*6 +  6*Nbre_thread_actif;
    E_Int* ipt_ind_CLgmres        = ipt_ind_CL         + (ithread-1)*6 + 12*Nbre_thread_actif;
    E_Int* ipt_shift_lu           = ipt_ind_CL         + (ithread-1)*6 + 18*Nbre_thread_actif;
    E_Int* ipt_ind_dm_socket      = ipt_ind_dm_omp     + (ithread-1)*12;
    E_Int* ipt_ind_dm_omp_thread  = ipt_ind_dm_socket  + 6;

    E_Int ipt_ind_sdm_thread[6];
    E_Int ipt_ind_coe_thread[6];
    E_Int ipt_ind_grad_thread[6];

    E_Int* ipt_nidom_loc, nb_subzone;

    /****************************************************
      -----Boucle sous-iteration
    ****************************************************/
    //cout << "nticfg= " << nitcfg << endl;
    if( nitcfg == 1){
      for (E_Int nd = 0; nd < nidom; nd++) { //mise a jour metric et vent ale zone cart et 3dhom(3dfull et 2d a la volee)
	if(param_int[nd][LALE]==1 && param_int[nd][ITYPZONE]!= 4 ) { //maillage indeformable structure
	  mjr_ale_3dhomocart_(nd, param_int[nd] ,   param_real[nd]   ,
			      socket             ,  Nbre_socket       , ithread_sock        , thread_parsock,
			      ipt_ind_dm_socket  , ipt_topology_socket,
			      iptx[nd]           , ipty[nd]           , iptz[nd]            ,
			      ipti[nd]           , iptj[nd]           , iptk[nd]            ,
			      ipti0[nd]          , iptj0[nd]          , iptk0[nd]           , iptvol[nd] ,
			      iptventi[nd]       , iptventj[nd]       , iptventk[nd]        );
	  //modifier mjr_ale_3dhomocart_ pour faire sauter barrier
#pragma omp barrier
	}
      }//zone

      // calcul metric si maillage deformable structure
#include           "../../FastS/FastS/Metric/cp_metric.cpp"
    }


    /// calcul de ndim

    E_Int ndim = 0;
    //---------------------------------------------------------------------
    // -----Boucle sur num.les domaines de la configuration
    // ---------------------------------------------------------------------
    E_Int shift_zone=0; E_Int shift_wig=0; E_Int shift_coe=0; E_Int nd_current=0;E_Int shift_rk4=0;  E_Int shift_grad=0; E_Int shift_mu=0; E_Int ncells=0;
    if  (param_int[0][EXPLOC] == 0 and param_int[0][RK] == 4 ){ //or param_int[0][EXPLOC] == 0 and param_int[0][RK] == 5)
      for (E_Int nd = 0; nd < nidom; nd++){ shift_rk4 = shift_rk4 + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ]; }
    }
    E_Float rhs_end=0;

    shift_rk4 = shift_rk4*(nitcfg - 1);
    //cout << "shift_rk4= " << shift_rk4 << endl;
	

    if (param_int[0][IMPLICITSOLVER] == 1 && layer_mode == 1) { ipt_norm_kry[ithread-1]=0.; }
    
    //calcul du sous domaine a traiter par le thread
    E_Int nbtask = ipt_omp[nitcfg-1]; 
    E_Int ptiter = ipt_omp[nssiter+ nitcfg-1];

    //printf("nbtaaask %d %d \n", nbtask, ptiter);
    for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
      E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
      E_Int nd = ipt_omp[ pttask ];
 
      E_Int cycl = param_int[nd][NSSITER]/param_int[nd][LEVEL];

      shift_zone=0; shift_wig=0; shift_coe=0; shift_grad=0; ndim =0;
      for (E_Int n = 0; n < nd; n++)
        {
          shift_zone = shift_zone + param_int[n][ NDIMDX ]*param_int[n][ NEQ ];
          shift_coe  = shift_coe  + param_int[n][ NDIMDX ]*param_int[n][ NEQ_COE ];
          if(param_int[n][ KFLUDOM ]==2){  shift_wig  = shift_wig  + param_int[n][ NDIMDX ]*3;}
          if(param_int[n][ITYPZONE ]==4){  shift_grad = shift_grad + param_int[n][ NDIMDX ]*param_int[n][ NEQ ]*3;}
          if (param_int[0][EXPLOC] == 1 and param_int[0][RK] == 3 and cycl != 4 and nitcfg%cycl==cycl/4){ ndim = ndim + param_int[nd][ NDIMDX ]*param_int[nd][ NEQ ]; }
        }
   
      E_Int lmin = 10;
      if (param_int[nd][ITYPCP] == 2) lmin = 4;

      if (param_int[nd][ITYPZONE] == 4){
#include       "../../FastP/FastP/Compute/rhs.cpp"
      }
      else
      {
	if (param_int[nd][IFLOW] == 4){
	  lmin =1;
#include       "../../FastLBM/FastLBM/Compute/cpp_files/lbm_collision_propagation.cpp"
	}
	else{
#include       "../../FastS/FastS/Compute/rhs.cpp"
	}
      }
      //E_Float fin_zone = omp_get_wtime();

      //if(ithread==1){cout <<"zone : "<< nd <<" "<< "temps= " << fin_zone - deb_zone <<" "<<"cycle =  " << cycl << endl;}
	   

    } //Fin boucle sur zones pour calcul RHS
          
#ifdef _WIN32
#pragma omp barrier
#endif
    //
    //timer pour omp "dynamique"
    //
#ifdef _OPENMP  
    rhs_end = omp_get_wtime();
#endif
    E_Int cpu_perthread = (ithread-1)*2 +(nitcfg-1)*Nbre_thread_actif*2;
    timer_omp[ cpu_perthread] += rhs_end - rhs_begin;
    //if(ithread==1) printf(" time rhs= %g %g  %d %d  \n",timer_omp[ cpu_perthread], rhs_end - rhs_begin , cpu_perthread, nitcfg);

    //#include   "../../FastS/FastS/Compute/lhs.cpp"
    //
    //finalisation timer pour omp "dynamique"
    //
#ifdef _OPENMP  
    E_Float     lhs_end = omp_get_wtime();
#else  
    E_Float     lhs_end = 0.;
#endif
    timer_omp[ cpu_perthread +1 ] += lhs_end- rhs_end;
    // if(ithread==1) printf(" time lhs= %g \n",lhs_end - rhs_end );
          

  } // Fin zone // omp
