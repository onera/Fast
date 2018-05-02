/*    
    Copyright 2013-2018 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

# include "fastS.h"
# include "param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

// 
// 
// Souszonelist (C layer)
// 
// 
void K_FASTS::distributeThreads_c( E_Int**& param_int, E_Int**& ipt_ind_dm, 
                                   E_Int& nidom  , E_Int& nssiter , E_Int& mx_omp_size_int  , E_Int& nstep, E_Int& nitrun, E_Int& display)
{
  // calcul  nombre souszone
  E_Int mxzone=0;
  for (E_Int nd = 0; nd < nidom; nd++)
     {  
       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
       E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   
       mxzone +=nb_subzone;
     } // loop zone

  FldArrayI tab_nijk(mxzone*3);   E_Int* ipt_nijk    = tab_nijk.begin();
  FldArrayI tab_ndimdx(mxzone);   E_Int* ipt_ndimdx  = tab_ndimdx.begin();
  FldArrayI tab_nozone(mxzone);   E_Int* ipt_nozone  = tab_nozone.begin();
  FldArrayI tab_nosszone(mxzone); E_Int* ipt_nosszone= tab_nosszone.begin();

  FldArrayI tab_size_subzone(nidom);   E_Int*  size_subzone = tab_size_subzone.begin();

  FldArrayI newtab_nijk(mxzone*3);   E_Int* ipt_nijk_new    = newtab_nijk.begin();
  FldArrayI newtab_ndimdx(mxzone);   E_Int* ipt_ndimdx_new  = newtab_ndimdx.begin();
  FldArrayI newtab_nozone(mxzone);   E_Int* ipt_nozone_new  = newtab_nozone.begin();
  FldArrayI newtab_nosszone(mxzone); E_Int* ipt_nosszone_new= newtab_nosszone.begin();

  E_Int c = 0;
  E_Int boom=0;
  E_Int check_size;
  E_Int ndimt=0;
  //calcul dimension et nombre souszone
  for (E_Int nd = 0; nd < nidom; nd++)
     {  
       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
       E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   

       size_subzone[nd] = 0;

       //Initialisation pointer iteration si zone skipper
       if(nb_subzone==0){ E_Int Ptomp =  param_int[nd][PT_OMP]; param_int[nd][Ptomp +nstep ] =  param_int[nd][Ptomp +nstep -1];}

       for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
         {  
            E_Int shift  = (nstep-1)*6*param_int[nd][MXSSDOM_LU]  + 6*nd_subzone;

            E_Int* nijk     = ipt_nijk     + 3*c;
            E_Int* ndimdx   = ipt_ndimdx   +   c;
            E_Int* nozone   = ipt_nozone   +   c;
            E_Int* nosszone = ipt_nosszone +   c;

            nijk[0]     = ipt_ind_dm[nd][1+shift];
            nijk[1]     = ipt_ind_dm[nd][3+shift];
            nijk[2]     = ipt_ind_dm[nd][5+shift];
            ndimdx[0]   = nijk[0]*nijk[1]*nijk[2];
            nozone[0]   = nd;
            nosszone[0] = nd_subzone;
            ndimt      += ndimdx[0];
            c +=1;
             
          }//souszone
     } // loop zone 

  E_Int Nbzones = c;

  //
  //
  //blindage pour lu_local si pas de travail a cette iteration
  //
  //
  if(ndimt==0)
  {
    for (E_Int c = 0; c < Nbzones; c++)
      {  
        E_Int* nijk         = ipt_nijk_new     + 3*c;
        E_Int* nozone_new   = ipt_nozone_new   +   c;
        E_Int* nosszone_new = ipt_nosszone_new +   c;

        E_Int  No_zone      = nozone_new[0];
        E_Int  No_sszone    = nosszone_new[0];
        E_Int* ndimdx       = ipt_ndimdx_new   +   c;

        E_Int* ipt_nidom_loc = ipt_ind_dm[No_zone] + param_int[No_zone][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
        E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   

        E_Int Ptomp =  param_int[No_zone][PT_OMP];

        if(nstep==1) 
        { //initilatisation pointer 1ere iter (les suivantes sont initilisees a lfin de cette routine)
                                          //    shift pointer iter | shift numa
          param_int[No_zone][Ptomp] =  Ptomp +       nssiter          + 1;
        } 
        //pointer pour iteration =nstep
        E_Int  PtrIterOmp    = param_int[No_zone][Ptomp +nstep -1];   

        //pointer pour la subzone courante a iteration =nstep
        E_Int  PtZoneomp    = PtrIterOmp + nb_subzone + size_subzone[No_zone];

        E_Int check = PtZoneomp - param_int[No_zone][PT_OMP  ];
        if(check > mx_omp_size_int && boom ==0) {boom =1; check_size =check;}
       
        param_int[No_zone][PtrIterOmp + No_sszone] = PtZoneomp;

        for (E_Int th = 0; th < __NUMTHREADS__; th++) { param_int[No_zone][ PtZoneomp + th ] = -2;} //Thread inactif
     }
   return;
  }


  //Tri par taille decroissante des zone
  for (E_Int c = 0; c < Nbzones; c++)
     {  

       E_Int ndimdx_max = 0; E_Int c_tg =-1;
       for (E_Int c1 = 0; c1 < Nbzones; c1++)
       {  
            E_Int* ndimdx = ipt_ndimdx  + c1;

            if(ndimdx[0] > ndimdx_max) { ndimdx_max= ndimdx[0]; c_tg = c1;}
       }// loop recherche plus grosse zone

       
       E_Int* nijk     = ipt_nijk     + 3*c_tg;
       E_Int* ndimdx   = ipt_ndimdx   +   c_tg;

       E_Int* nozone   = ipt_nozone   +   c_tg;
       E_Int* nosszone = ipt_nosszone +   c_tg;

       //printf("zone %d %d %d %d %d \n", c, ndimdx_max, c_tg, nozone[0],  nosszone[0]);

       E_Int* nijk_new     = ipt_nijk_new     + 3*c;
       E_Int* ndimdx_new   = ipt_ndimdx_new   +   c;
       E_Int* nozone_new   = ipt_nozone_new   +   c;
       E_Int* nosszone_new = ipt_nosszone_new +   c;

       nijk_new[0]     = nijk[0];
       nijk_new[1]     = nijk[1];
       nijk_new[2]     = nijk[2];
       ndimdx_new[0]   = ndimdx[0];
       nozone_new[0]   = nozone[0];
       nosszone_new[0] = nosszone[0];
       ndimdx[0]       = -1;
     } // loop zone


  E_Int cells_tg = 0;
  E_Int flui_tg  = 0;
  E_Int fluj_tg  = 0;
  E_Int fluk_tg  = 0;
  
  FldArrayI inddm_zone(Nbzones*6);   E_Int* ipt_inddm_zone    = inddm_zone.begin();

  for (E_Int c = 0; c < Nbzones; c++)
     {  
        E_Int* nijk   = ipt_nijk_new     + 3*c;
        E_Int* ind_dm = ipt_inddm_zone   + 6*c;

        ind_dm[0]= 1; ind_dm[2]= 1; ind_dm[4]= 1;
        ind_dm[1]= nijk[0]; ind_dm[3]= nijk[1]; ind_dm[5]= nijk[2];

	E_Int size_c = nijk[0]*nijk[1]*nijk[2];
	E_Int size_i =(nijk[0]+1)*nijk[1]*nijk[2];
	E_Int size_j =(nijk[1]+1)*nijk[0]*nijk[2];
	E_Int size_k = 0;
        if(nijk[2] != 1) size_k = (nijk[2]+1)*nijk[0]*nijk[1];

        cells_tg += size_c;
        flui_tg  += size_i;
        fluj_tg  += size_j;
        fluk_tg  += size_k;
     }  

    cells_tg /= __NUMTHREADS__;
    flui_tg  /= __NUMTHREADS__;
    fluj_tg  /= __NUMTHREADS__;
    fluk_tg  /= __NUMTHREADS__;

  // creation et initialisation remainder
  FldArrayI tab_remaind(__NUMTHREADS__);   E_Int* remaind  = tab_remaind.begin();

  FldArrayI tab_topo_lu(3);   E_Int* topo_lu  = tab_topo_lu.begin();
  FldArrayI tab_ind_dm_th(6); E_Int* ind_dm_th= tab_ind_dm_th.begin();

  FldArrayI tab_dim_i(__NUMTHREADS__);   E_Int* dim_i  = tab_dim_i.begin();
  FldArrayI tab_dim_j(__NUMTHREADS__);   E_Int* dim_j  = tab_dim_j.begin();
  FldArrayI tab_dim_k(__NUMTHREADS__);   E_Int* dim_k  = tab_dim_k.begin();

  FldArrayI tab_numa_socket(NBR_SOCKET);   E_Int* numa_socket  = tab_numa_socket.begin();

  for (E_Int th = 0; th < __NUMTHREADS__ ; th++){ remaind[th] = cells_tg;}

  //
  E_Int cells_tg_save = cells_tg;
  E_Float      marge  = 1.005;

  for (E_Int c = 0; c < Nbzones; c++)
     {  
        E_Int* ind_dm       = ipt_inddm_zone   + 6*c;
        E_Int* nijk         = ipt_nijk_new     + 3*c;
        E_Int* nozone_new   = ipt_nozone_new   +   c;
        E_Int* nosszone_new = ipt_nosszone_new +   c;

        E_Int  No_zone      = nozone_new[0];
        E_Int  No_sszone    = nosszone_new[0];
        E_Int* ndimdx       = ipt_ndimdx_new   +   c;

        E_Int* ipt_nidom_loc = ipt_ind_dm[No_zone] + param_int[No_zone][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
        E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   

        //printf("c Nozone %d %d %d %d %d %d %d \n", c, No_zone, ind_dm[0], ind_dm[1], ind_dm[2], ind_dm[3], ndimdx[0] );

        E_Int Ptomp =  param_int[No_zone][PT_OMP];

        if(nstep==1) 
        { //initilatisation pointer 1ere iter (les suivantes sont initilisees a lfin de cette routine)
                                          //    shift pointer iter | shift numa
          param_int[No_zone][Ptomp] =  Ptomp +       nssiter          + 1;
          //initilisation compteur Numa
           for (E_Int socket = 0; socket < NBR_SOCKET; socket++) { numa_socket[socket] = 0;} 
        } 
        //pointer pour iteration =nstep
        E_Int  PtrIterOmp    = param_int[No_zone][Ptomp +nstep -1];   

        //pointer pour la subzone courante a iteration =nstep
        E_Int  PtZoneomp    = PtrIterOmp + nb_subzone + size_subzone[No_zone];

        //printf("Ptomp= %d , ptiter= %d , pt_subzone= %d , nstep= %d \n ",Ptomp, PtrIterOmp, PtZoneomp ,   nstep);
        E_Int check = PtZoneomp - param_int[No_zone][PT_OMP  ];
        if(check > mx_omp_size_int && boom ==0) {boom =1; check_size =check;}
       
        param_int[No_zone][PtrIterOmp + No_sszone] = PtZoneomp;

	E_Int size_c = nijk[0]*nijk[1]*nijk[2];
	E_Int size_i =(nijk[0]+1)*nijk[1]*nijk[2];
	E_Int size_j =(nijk[1]+1)*nijk[0]*nijk[2];
	E_Int size_k = 0;
        if(nijk[2] != 1) size_k = (nijk[2]+1)*nijk[0]*nijk[1];

        E_Int Nthreads;
        //printf("so %d %f%d \n ", size_c, cells_tg*marge, c);

        if(size_c <= cells_tg*marge)
        {
	    E_Int th_current = 0;
            E_Int lgo = 1;
	    while(lgo==1)
              { //print 'th_current',th_current,remaind[str(th_current)],size_c/marge, size_c

                if(size_c/marge <=  remaind[th_current]) lgo = 0;
                if(lgo==1 && th_current == __NUMTHREADS__-1)
                   {
                     lgo = 0;
                     //on cherche le thread le moins charge pour attribuer le residu
                     E_Int th_max =-100000;
                     for (E_Int th = 0; th < __NUMTHREADS__; th++)
                       {  
                        if(remaind[th] > th_max){ th_max = remaind[th]; th_current = th;}
                        //print 'traitememnt residu affecte au thread', th_current, remaind[str(th)], th
                       }
                   }
                if(lgo==1) th_current +=1;
              }// while go
  
            remaind[th_current] =  remaind[th_current] - size_c; 
             
            //Carte threads actifs       
            for (E_Int th = 0; th < __NUMTHREADS__; th++) { param_int[No_zone][ PtZoneomp + th ] = -2;} //Thread inactif

            param_int[No_zone][ PtZoneomp + th_current                 ] =  0; //le Thread "0" (pas au sens omp: indirection entre omp_numthread et ithread)) prend le job
            //Nb threads actifs 
            Nthreads = 1;
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__   ] =  Nthreads;
            //topology omp
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +1] =  1;
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +2] =  1;
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +3] =  1;
            //indice sous-domaine omp

            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +4 ] =  ind_dm[0];
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +5 ] =  ind_dm[1];
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +6 ] =  ind_dm[2];
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +7 ] =  ind_dm[3];
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +8 ] =  ind_dm[4];
            param_int[No_zone][ PtZoneomp + __NUMTHREADS__ +9 ] =  ind_dm[5];

            E_Int check = PtZoneomp + __NUMTHREADS__ +9 - param_int[No_zone][PT_OMP  ];
            if(check > mx_omp_size_int  && boom ==0) {boom =1; check_size =check;}


           if ((display==1 && nitrun ==1 && nstep==1) || display ==2)
           {
             printf("souszone %d de la zone %d calcule par %d  Threads,  nstep= %d \n", c, No_zone, Nthreads, nstep);
             printf("-- Thread= %d : fenetre=( %d, %d, %d, %d, %d, %d) \n", th_current, ind_dm[0],ind_dm[1],ind_dm[2],ind_dm[3],ind_dm[4],ind_dm[5]);
           }
        }
        else
        {
            E_Int cells_tg = cells_tg_save;
            Nthreads =  size_c/cells_tg;
            if(size_c%cells_tg != 0) Nthreads +=1;
            Nthreads = K_FUNC::E_min( Nthreads,  __NUMTHREADS__ );

            //verif place dispo sur les dernier threads
            
            cells_tg     = 0;
            E_Int cells_tg_min =10000000; E_Int th_min;

            for (E_Int th_check = __NUMTHREADS__-1; th_check > __NUMTHREADS__ -Nthreads-1; th_check--)
               { cells_tg  +=  remaind[th_check];
                 if(remaind[th_check] <  cells_tg_min)
                    { 
                      cells_tg_min = remaind[th_check];
                      th_min = th_check;
                    }
               }
            cells_tg  /=  Nthreads;

            if(cells_tg == cells_tg_min) th_min =-1;
            cells_tg  =  cells_tg_min*marge;

            E_Int lmin;
            if(param_int[No_zone][ITYPCP] == 2 ) lmin =  4;
            else  lmin = 10;

            FldArrayI tab_dim_ijk(3*Nthreads);   E_Int* dim_ijk  = tab_dim_ijk.begin();
            FldArrayI tab_res_ijk(3*Nthreads);   E_Int* res_ijk  = tab_res_ijk.begin();


            E_Int search = 1;
            while(search == 1)
              {
                //#verif place dispo sur les dernier threads
                //#for th_check in range(OMP_NUM_THREADS-1,OMP_NUM_THREADS-Nthreads-1,-1):
                //#   print 'remain', remaind[str(th_check)], th_check

                search = 0;

                //#print 'ind_zone', inddm_zones[c], Nthreads
                E_Int ithread = 1;
                indice_boucle_lu_(c, ithread, Nthreads, lmin, ind_dm, topo_lu, ind_dm_th );
                //#print "topo_LLLUUU", topo_lu

                
                for (E_Int dir = 0; dir < 3; dir++){
                  for (E_Int l = 0; l < topo_lu[dir]; l++){

                        if     (dir==0){dim_i[l] =  nijk[dir]/topo_lu[dir]; if(l < nijk[dir]%topo_lu[dir]){dim_i[l] += 1;}  }
                        else if(dir==1){dim_j[l] =  nijk[dir]/topo_lu[dir]; if(l < nijk[dir]%topo_lu[dir]){dim_j[l] += 1;}  }
                        else           {dim_k[l] =  nijk[dir]/topo_lu[dir]; if(l < nijk[dir]%topo_lu[dir]){dim_k[l] += 1;}  }
                    }
                  }

                E_Int res = dim_i[0]*dim_j[0]*dim_k[0] - cells_tg;
                E_Float sign = 1.;
                if(res < 0) sign = -1.; //block trop petit % la cible

                //#print 'size bloc=',dims_i[0]*dim_j[0]*dim_k[0 , 'sizetg=', cells_tg
            
                //sous bloc trop grand
                if(res*sign > 0)
                  {
                    E_Int compteur =0;
                    E_Int go = 0;
                    while(res*sign > 0. && go==0 && compteur < 6)
                       {
                        go      = 1;
                        E_Int deriv_i= K_FUNC::E_min(1, (topo_lu[0]-1) );
                        E_Int deriv_j= K_FUNC::E_min(1, (topo_lu[1]-1) );
                        E_Int deriv_k= K_FUNC::E_min(1, (topo_lu[2]-1) );

                        E_Int cout_i = dim_j[0]*dim_k[0]*deriv_i;
                        E_Int cout_j = dim_i[0]*dim_k[0]*deriv_j;
                        E_Int cout_k = dim_j[0]*dim_i[0]*deriv_k;
                        if( cout_i+cout_j+cout_k < res*sign)
                          {
                            dim_i[0] -=deriv_i*sign; 
                            dim_j[0] -=deriv_j*sign;
                            dim_k[0] -=deriv_k*sign;
                            res      -= (cout_i+cout_j+cout_k)*sign;
                            if(cout_i+cout_j+cout_k != 0) go = 0;
                          }
                        else if( cout_j+cout_k < res*sign)
                          {
                            dim_j[0] -=deriv_j*sign;
                            dim_k[0] -=deriv_k*sign;
                                 res -= (cout_j+cout_k)*sign;
                            if(cout_j+cout_k != 0) go = 0;
                          }
                        else if(cout_j+cout_i < res*sign)
                          {
                            dim_j[0] -=deriv_j*sign; 
                            dim_i[0] -=deriv_i*sign;
                            res      -= (cout_j+cout_i)*sign;
                            if(cout_j+cout_i != 0) go = 0;
                          }
                        else if(cout_k+cout_i < res*sign)
                          {
                            dim_k[0] -=deriv_k*sign; 
                            dim_i[0] -=deriv_i*sign;
                            res      -= (cout_k+cout_i)*sign;
                            if( cout_k+cout_i != 0) go = 0;
                          }
                        else if(cout_k < res*sign)
                          {
                            dim_k[0] -=deriv_k*sign;
                            res      -= cout_k*sign;
                            if(cout_k != 0) go = 0;
                          }
                        else if(cout_j < res*sign)
                          {
                            dim_j[0] -=deriv_j*sign;
                            res      -= cout_j*sign;
                            if(cout_j != 0) go = 0;
                          }
                        else if(cout_i < res*sign)
                          {
                            dim_i[0] -=deriv_i*sign;
                            res      -= cout_i*sign;
                            if(cout_i != 0) go = 0;
                          }
                        compteur +=1;
                       }//while go
                  } // if res

              for (E_Int dir = 0; dir < 3; dir++){
                   if(topo_lu[dir] != 1)
                   {
                      //on repartie la place restante entre les thread restant
                      if(th_min != -1)
                      {
                          E_Int test;
                          if     (dir==0) test = dim_i[0];
                          else if(dir==1) test = dim_j[0];
                          else            test = dim_k[0];

                          E_Int dim_loc = (nijk[dir]-test)/(topo_lu[dir]-1);
                          E_Int res_loc = (nijk[dir]-test)%(topo_lu[dir]-1);
                          for( E_Int l=1; l < topo_lu[dir]; l++){
                              if     (dir==0){dim_i[l] =  dim_loc; if(l < res_loc){dim_i[l] += 1;}  }
                              else if(dir==1){dim_j[l] =  dim_loc; if(l < res_loc){dim_j[l] += 1;}  }
                              else           {dim_k[l] =  dim_loc; if(l < res_loc){dim_k[l] += 1;}  }
                            }
                      }
                      //#on affecte la taille du premier bloc thraed a tous les autres bloc
                      else
                      {
                          for( E_Int l=1; l < topo_lu[dir]; l++){ 
                              if     (dir==0)dim_i[l] =  dim_i[0]; 
                              else if(dir==1)dim_j[l] =  dim_j[0];
                              else           dim_k[l] =  dim_k[0];
                            }
                          //on determine la taille du dernier bloc 
                          if     (dir==0) dim_i[ topo_lu[dir]-1 ] = nijk[dir] - (topo_lu[dir]-1)*dim_i[0];
                          else if(dir==1) dim_j[ topo_lu[dir]-1 ] = nijk[dir] - (topo_lu[dir]-1)*dim_j[0];
                          else            dim_k[ topo_lu[dir]-1 ] = nijk[dir] - (topo_lu[dir]-1)*dim_k[0];
                      }
                   }
                 } // loop dir

              for (E_Int dir = 0; dir < 3; dir++){
                for( E_Int l=0; l < topo_lu[dir]; l++){
                      E_Int test;
                      if     (dir==0) test = dim_i[l];
                      else if(dir==1) test = dim_j[l];
                      else            test = dim_k[l];
                      if(search==0 && test < lmin &&  nijk[dir] != test)
                      {  
                        search = 1;
                        if     (dir ==0 ) Nthreads -=  topo_lu[1]*topo_lu[2];
                        else if(dir ==1 ) Nthreads -=  topo_lu[0]*topo_lu[2];
                        else if(dir ==2 ) Nthreads -=  topo_lu[0]*topo_lu[1];
                      }
                    }//loop topo
                 } // loop dir
              }// while search

            //Carte threads actifs       
            for (E_Int th = 0; th < __NUMTHREADS__; th++) { param_int[No_zone][ PtZoneomp + th ] = -2;} //Thread inactif

            //Recherche thread disponible
            E_Int th = 0;
            for (E_Int k = 0; k < topo_lu[2]; k++){
              for (E_Int j = 0; j < topo_lu[1]; j++){
                for (E_Int i = 0; i < topo_lu[0]; i++){
		  E_Int th_current = 0;
                  E_Int size_th    = dim_i[i]*dim_j[j]*dim_k[k];
                  E_Int size_loc   = size_th;

	          while(size_loc/marge > remaind[th_current] && th_current < __NUMTHREADS__-1) th_current +=1;
                
                  if(th_current == __NUMTHREADS__-1)
                  {
                     E_Int blind = -100000;
                     for (E_Int l = 0; l < __NUMTHREADS__; l++){
                        if(remaind[l] > blind)
                        { blind = remaind[l];
                          th_current = l;
                        }
                     }
                  }
                  //tentative optim distrib
                  if(th_min != -1) { th_current = th_min; th_min = -1; }

                  //on estime a quelle socket Numa la zone appartient  ( a blinder...)
                  E_Int thread_parsoc = __NUMTHREADS__/NBR_SOCKET;
                  if (thread_parsoc==0) thread_parsoc=1;
                  E_Int socket = th_current/thread_parsoc;
                  if (socket >= NBR_SOCKET) socket=0;
                  numa_socket[socket] += size_loc;

                  //printf("th_current %d %d %d %d %d %d \n", th_current, remaind[th_current], i,j,k, nstep  );

                  remaind[th_current]                          =  remaind[th_current] - size_loc;
                  param_int[No_zone][ PtZoneomp + th_current ] = th;                                    //Thread actif
                  th +=1;
                  //print 'thread=', th, size_loc, remaind[str(th_current)], th_current
                }//loop i
              }//loop j
            }//loop k

            if(nstep ==1)
            {E_Int numa_max=0; E_Int sock_tg =0;
             for (E_Int socket = 0; socket < NBR_SOCKET; socket++) { if( numa_socket[socket] > numa_max) { numa_max= numa_socket[socket]; sock_tg= socket;} } 
             param_int[No_zone][Ptomp + nssiter]= sock_tg;
            }
            //Nb threads actifs
            param_int[No_zone][ PtZoneomp +  __NUMTHREADS__   ] =  topo_lu[0]*topo_lu[1]*topo_lu[2];
            //topology omp
            param_int[No_zone][ PtZoneomp +  __NUMTHREADS__ +1] =  topo_lu[0];
            param_int[No_zone][ PtZoneomp +  __NUMTHREADS__ +2] =  topo_lu[1];
            param_int[No_zone][ PtZoneomp +  __NUMTHREADS__ +3] =  topo_lu[2];

            th = 0;
            E_Int kstart = ind_dm[4];
            for (E_Int k = 0; k < topo_lu[2]; k++){
              E_Int jstart = ind_dm[2];
              for (E_Int j = 0; j < topo_lu[1]; j++){
                E_Int istart = ind_dm[0];
                for (E_Int i = 0; i < topo_lu[0]; i++){


                    ind_dm_th[0] = istart;
                    ind_dm_th[1] = istart + dim_i[i] - 1;
                    istart       = ind_dm_th[1] + 1;

                    ind_dm_th[2] = jstart;
                    ind_dm_th[3] = jstart + dim_j[j] - 1;

                    ind_dm_th[4] = kstart;
                    ind_dm_th[5] = kstart + dim_k[k] - 1;
                   

                    param_int[No_zone][PtZoneomp + __NUMTHREADS__ +4 + th*6] = ind_dm_th[0];
                    param_int[No_zone][PtZoneomp + __NUMTHREADS__ +5 + th*6] = ind_dm_th[1];
                    param_int[No_zone][PtZoneomp + __NUMTHREADS__ +6 + th*6] = ind_dm_th[2];
                    param_int[No_zone][PtZoneomp + __NUMTHREADS__ +7 + th*6] = ind_dm_th[3];
                    param_int[No_zone][PtZoneomp + __NUMTHREADS__ +8 + th*6] = ind_dm_th[4];
                    param_int[No_zone][PtZoneomp + __NUMTHREADS__ +9 + th*6] = ind_dm_th[5];

                    th     +=1;
                }//loop i
               jstart  =   ind_dm_th[3]+ 1;
              }//loop j
             kstart  = ind_dm_th[5] + 1;
            }//loop k

            E_Int check = PtZoneomp + __NUMTHREADS__ +9 +th*6 - param_int[No_zone][PT_OMP  ];
            if(check > mx_omp_size_int  && boom ==0 ) {boom =1; check_size =check;}

           if ((display==1 && nitrun ==1 && nstep==1) || display==2)
           { 
             E_Int Nthreads = param_int[No_zone][ PtZoneomp + __NUMTHREADS__ ];

             printf("souszone %d de la zone %d calcule par %d  Threads,  nstep= %d \n", c, No_zone, Nthreads, nstep);

             E_Int th=0;
             for (E_Int th_current = 0; th_current < __NUMTHREADS__; th_current++){

                 if( param_int[No_zone][ PtZoneomp + th_current ] != -2)
                 {
                   E_Int* ind_dm = param_int[No_zone] + PtZoneomp + __NUMTHREADS__ +4 + th*6;

                   printf("-- Thread= %d : fenetre=( %d, %d, %d, %d, %d, %d) \n", th_current,
                                                 ind_dm[0],ind_dm[1],ind_dm[2],ind_dm[3],ind_dm[4],ind_dm[5]);
                   th +=1;
                 }
             }//loop
           }//if display

        }// if mono ou subzone

     size_subzone[No_zone] += 4 + __NUMTHREADS__ + 6*Nthreads;

       if(boom == 1 )
               {
                printf("------\n");
                printf("Error msg\n");
                printf("------\n");
                printf("resize MX_OMP_SIZE_INT. Present value= %d \n ", mx_omp_size_int);
                printf("Value must be at least larger than : %d \n ", check_size);
                printf("Just after the modules import of userscript.py, add the following python command:\n");
                printf("#\n");
                printf("#\n");
                printf("Fast.FastI.MX_OMP_SIZE_INT= %d\n ", check_size+20);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }

     }//loop zone

   //mise a jour Pointer pour stockage nstep+1
   if(nstep != nssiter)
   {
     for (E_Int c = 0; c < Nbzones; c++)
     {  
        E_Int* nozone_new   = ipt_nozone_new   +   c;
        E_Int  No_zone      = nozone_new[0];

        E_Int* ipt_nidom_loc = ipt_ind_dm[No_zone] + param_int[No_zone][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
        E_Int nb_subzone     = ipt_nidom_loc [nstep-1];


        E_Int Ptomp =  param_int[No_zone][PT_OMP];
        param_int[No_zone][Ptomp +nstep ] =  param_int[No_zone][Ptomp +nstep -1] + nb_subzone + size_subzone[No_zone];

        E_Int check_size = param_int[No_zone][Ptomp +nstep ] - param_int[No_zone][PT_OMP  ];
        if(check_size > mx_omp_size_int)
               {
                printf("------\n");
                printf("Error msg\n");
                printf("------\n");
                printf("resize MX_OMP_SIZE_INT. Present value= %d \n ", mx_omp_size_int);
                printf("Value must be at least larger than : %d \n ", check_size);
                printf("Just after the modules import of userscript.py, add the following python command:\n");
                printf("#\n");
                printf("#\n");
                printf("Fast.FastI.MX_OMP_SIZE_INT= %d\n ", check_size+20);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }
     }  
   }

  if ( (display==1 && nitrun ==1 && nstep == 1) || display ==2)
  {  printf("OMP balancing: positive values --> slow thread. Target cells per Thread: %d Nitrun= %d, nstep= %d \n", cells_tg_save, nitrun , nstep);
     for (E_Int th_current = 0; th_current < __NUMTHREADS__; th_current++){
        printf("desequilibre thread, %d, : %f percent \n", th_current, (float(-remaind[th_current])/cells_tg)*100.);
     }
  }
 
  return;
}


// -----------------------------------------------------------------------------------
// 
// distributeThreads (python layer)
// 
// -----------------------------------------------------------------------------------
PyObject* K_FASTS::distributeThreads(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; 
  E_Int nssiter; E_Int nstep; E_Int nitrun; E_Int mx_omp_size_int; E_Int display;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOlllll", &zones , &metrics, &nstep, &nssiter, &nitrun, &mx_omp_size_int, &display)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOiiiii", &zones , &metrics, &nstep, &nssiter, &nitrun, &mx_omp_size_int, &display)) return NULL; 
#endif
  
  E_Int nidom = PyList_Size(zones);

  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; 
  ipt_param_int     = new E_Int*[nidom*2];
  ipt_ind_dm        = ipt_param_int   + nidom;

  vector<PyArrayObject*> hook;

  for (E_Int nd = 0; nd < nidom; nd++)
  {    
    PyObject* zone = PyList_GetItem(zones, nd);
       
    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);

    // get metric
    PyObject* metric     = PyList_GetItem(metrics, nd); // metric du domaine i
    ipt_ind_dm[nd]       = K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );
  }

  distributeThreads_c( ipt_param_int , ipt_ind_dm, nidom, nssiter, mx_omp_size_int, nstep, nitrun, display );

  delete [] ipt_param_int;

  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}
