/*    
    Copyright 2013-2023 Onera.

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

# include "FastC/fastc.h"
# include "FastC/param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

// 
// 
// Souszonelist (C layer)
// 
// 
void K_FASTC::distributeThreads_c( E_Int**& param_int, E_Float**& param_real, E_Int**& ipt_ind_dm,
                                   E_Int& nidom  , E_Int* ipt_dtloc , E_Int& mx_omp_size_int  , E_Int& nstep, E_Int& nitrun, E_Int& display )
{
#ifdef NB_SOCKET
  E_Int nb_socket=NB_SOCKET;
#else
  E_Int nb_socket=1;
#endif
#ifdef CORE_PER_SOCK
  E_Int core_per_socket=CORE_PER_SOCK;
#else
  E_Int core_per_socket=__NUMTHREADS__;
#endif

  E_Int nssiter  = ipt_dtloc[0];
  E_Int omp_mode = ipt_dtloc[8];
  E_Int shift_omp= ipt_dtloc[11];
  E_Int* ipt_omp = ipt_dtloc + shift_omp;

  E_Int lmin;

  if(omp_mode==0)
  {
    //calcul nbe sszone total
    E_Int mxzone=0;
    for (E_Int nd = 0; nd < nidom; nd++)
       {
          E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
          E_Int nb_subzone    = ipt_nidom_loc [nstep-1]; 
          mxzone += nb_subzone;
       }

    if(nstep==1){ipt_omp[ nssiter ]= 2*nssiter;}

    ipt_omp[ nstep-1 ]= mxzone;
    //printf("mxzone %d %d %d %d  \n", mxzone, 9 +  nssiter + nstep-1, nssiter, nstep);

    E_Int PtiterTask = ipt_omp[nssiter + nstep-1];

    E_Int c=0;
    for (E_Int nd = 0; nd < nidom; nd++)
       {
          E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
          E_Int nb_subzone    = ipt_nidom_loc [nstep-1]; 
           

#include "FastC/HPC_LAYER/SIZE_MIN.cpp"

          for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
          {
            E_Int* ipt_ind_dm_loc   = ipt_ind_dm[nd]  + (nstep-1)*6*param_int[nd][ MXSSDOM_LU ] + 6*nd_subzone;

            E_Int PtTask = PtiterTask +c*(6+ __NUMTHREADS__*7);
            ipt_omp[PtTask                     ] = nd;
            ipt_omp[PtTask +1                  ] = nd_subzone;
            ipt_omp[PtTask + __NUMTHREADS__ +2 ] = __NUMTHREADS__;

            E_Int* ipt_topo   = ipt_omp + PtTask + __NUMTHREADS__ +3;

            for (E_Int th = 1; th < __NUMTHREADS__ +1; th++)
            {
              E_Int* ind_dm = ipt_omp + PtTask + 6 + __NUMTHREADS__ +6*(th-1);

              E_Int Nbre_thread_actif_loc = __NUMTHREADS__;
              indice_boucle_lu_(nd, th, Nbre_thread_actif_loc , lmin, ipt_ind_dm_loc, ipt_topo, ind_dm);

              ipt_omp[PtTask + 1 + th  ] = th-1;

              //if( ind_dm[1] < ind_dm[0])
              // {
              //  ipt_omp[PtTask + 1 + th            ] = -2;
              //  ipt_omp[PtTask + __NUMTHREADS__ +2 ] = 1; //zone trop petite: 1 seul thread actif
              // }
              if ((display==1 && nitrun ==1 && nstep==1) || display ==2)
               {
                if(th==1){printf("souszone %d de la zone %d calcule par %d  Threads,  nstep= %d, cycle = %d \n", c, nd,ipt_omp[PtTask+ __NUMTHREADS__+2], nstep, param_int[nd][NSSITER]/param_int[nd][LEVEL] );}
                printf("-- Thread= %d : fenetre=( %d, %d, %d, %d, %d, %d) \n", th, ind_dm[0],ind_dm[1],ind_dm[2],ind_dm[3],ind_dm[4],ind_dm[5]);
               }
            }
          c+=1;
          } //loop sszones
       } //loop zones

   if(nstep != nssiter) 
      { 
       ipt_omp[ nssiter +nstep  ]= ipt_omp[ nssiter +nstep-1  ] + mxzone*(6+7*__NUMTHREADS__);

        E_Int check_size = ipt_omp[ nssiter +nstep  ] - shift_omp;
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
                printf("Fast.FastC.MX_OMP_SIZE_INT= %d\n ", check_size+20);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }
      }

   return;
  }//omp_mode0

  if( __NUMTHREADS__ <= core_per_socket)  nb_socket=1;

  //printf(" socket core %d %d \n", nb_socket, core_per_socket);
  
  /// obsolete (Ivan Mary)??
  FldArrayI tab_activ_core_per_socket(nb_socket);   E_Int* activ_core_per_socket = tab_activ_core_per_socket.begin();

  E_Int count = __NUMTHREADS__;
  for (E_Int s = 0; s < nb_socket; s++)
     {   if (count >= core_per_socket) {activ_core_per_socket[s]=  core_per_socket; count -= core_per_socket;}
         else activ_core_per_socket[s] = count;
     }
  /// obsolete??
          

 //// Determination du plus grand niveau en tps pour dtloc////
 //
 //
 E_Int level_max;
 if (param_int[0][EXPLOC]==0){ level_max = 0; } // Tous les schemas sauf dtlocal instationnaire
 else { E_Int denominateur_max = param_int[0][NSSITER]/4; level_max = log(denominateur_max)/log(2); } // dtlocal instationnaire


 //
 //
 //
 //// Boucle sur les niveaux en temps pour la distrib ////
 //
 //
 //
 for (E_Int lev = 0; lev < level_max + 1; lev ++)
 {
   // calcul  nombre souszone pour le niveau en temps lev
   E_Int zones_dtloc[nidom];
   E_Int nb_zones_dtloc =0;
   if (param_int[0][EXPLOC]==0)
     {  for (E_Int nd = 0; nd < nidom; nd ++) { zones_dtloc[ nb_zones_dtloc ] = nd; nb_zones_dtloc+= 1;} }
   else
     {
      for (E_Int nd = 0; nd < nidom; nd ++)
        { if(param_int[nd][LEVEL] == pow(2,lev) ) { zones_dtloc[ nb_zones_dtloc ] = nd; nb_zones_dtloc+= 1; } }
     }

   E_Int mxzone=0;
   for (E_Int i = 0; i < nb_zones_dtloc; i++)
    {
      E_Int nd = zones_dtloc[i];
      E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
      E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   
      mxzone += nb_subzone;
    } // loop zone

 
   if(nstep==1){ipt_omp[ nssiter ]= 2*nssiter;}

   FldArrayI tab_nijk(mxzone*3*2); E_Int* ipt_nijk    = tab_nijk.begin();
   FldArrayI tab_ndimdx(mxzone);   E_Int* ipt_ndimdx  = tab_ndimdx.begin();
   FldArrayI tab_nozone(mxzone);   E_Int* ipt_nozone  = tab_nozone.begin();
   FldArrayI tab_nosszone(mxzone); E_Int* ipt_nosszone= tab_nosszone.begin();

   //FldArrayI tab_size_task(mxzone);     E_Int*  size_task    = tab_size_task.begin();

   FldArrayI newtab_nijk(mxzone*3*2); E_Int* ipt_nijk_new    = newtab_nijk.begin();
   FldArrayI newtab_ndimdx(mxzone);   E_Int* ipt_ndimdx_new  = newtab_ndimdx.begin();
   FldArrayI newtab_nozone(mxzone);   E_Int* ipt_nozone_new  = newtab_nozone.begin();
   FldArrayI newtab_nosszone(mxzone); E_Int* ipt_nosszone_new= newtab_nosszone.begin();

   FldArrayF    tab_HPC_CUPS(mxzone); E_Float* ipt_HPC_CUPS    = tab_HPC_CUPS.begin();
   FldArrayF newtab_HPC_CUPS(mxzone); E_Float* ipt_HPC_CUPS_new= newtab_HPC_CUPS.begin();

   E_Int c = 0;
   E_Int ndimt=0;
   E_Float CupsMax =0;
   //calcul dimension et nombre souszone
   for (E_Int i = 0; i < nb_zones_dtloc ; i++)
     {  
       E_Int nd = zones_dtloc[i];
       E_Int* ipt_nidom_loc = ipt_ind_dm[nd] + param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
       E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   
       //cout << "zone " << nd << " " << nb_subzone << endl; 

       for (E_Int nd_subzone = 0; nd_subzone < nb_subzone; nd_subzone++)
         {  
            E_Int shift  = (nstep-1)*6*param_int[nd][MXSSDOM_LU]  + 6*nd_subzone;

            E_Int* nijk     = ipt_nijk     + 3*c;
            E_Int* ijk_start= ipt_nijk     + 3*c + 3*mxzone;
            E_Int* ndimdx   = ipt_ndimdx   +   c;
            E_Int* nozone   = ipt_nozone   +   c;
            E_Int* nosszone = ipt_nosszone +   c;

            ipt_HPC_CUPS[c] = param_real[nd][HPC_CUPS];
            if (ipt_HPC_CUPS[c] > CupsMax) CupsMax =  ipt_HPC_CUPS[c];


            ijk_start[0]= ipt_ind_dm[nd][0+shift];
            ijk_start[1]= ipt_ind_dm[nd][2+shift];
            ijk_start[2]= ipt_ind_dm[nd][4+shift];
            nijk[0]     = ipt_ind_dm[nd][1+shift]-ipt_ind_dm[nd][0+shift]+1;
            nijk[1]     = ipt_ind_dm[nd][3+shift]-ipt_ind_dm[nd][2+shift]+1;
            nijk[2]     = ipt_ind_dm[nd][5+shift]-ipt_ind_dm[nd][4+shift]+1;

            ndimdx[0]   = nijk[0]*nijk[1]*nijk[2];
            nozone[0]   = nd;
            nosszone[0] = nd_subzone;
            ndimt      += ndimdx[0];
            c +=1;
          }//souszone
     } // loop zone 

   //nombre global  de souszone a traiter a l'iter nstep
   ipt_omp[ nstep-1 ]= mxzone;

   E_Int boom       = 0;
   E_Int check_size = 0;
   E_Int Nbzones    = c;
   //blindage pour lu_local si pas de travail a cette iteration
   //
   //
   if(ndimt==0)
   {
    //ptitertask pointe sur l'iter precedente
    ipt_omp[ nssiter +nstep  ]= ipt_omp[ nssiter +nstep-1  ];
    
    continue; //on skip ce level
   }

   //Tri par taille decroissante des zone
   E_Float cout= 0.;
   for (E_Int c = 0; c < Nbzones; c++)
     {  

       E_Float ndimdx_max = 0; E_Int c_tg =-1;
       for (E_Int c1 = 0; c1 < Nbzones; c1++)
       {  
            E_Int* ndimdx = ipt_ndimdx  + c1;

            E_Float* hpc_cups = ipt_HPC_CUPS+   c1;

            //if(ndimdx[0] > ndimdx_max) { ndimdx_max= ndimdx[0]; c_tg = c1;}
            if( float(ndimdx[0])/hpc_cups[0] > ndimdx_max) { ndimdx_max= float(ndimdx[0])/hpc_cups[0]; c_tg = c1;}
       }// loop recherche plus grosse zone

       
       E_Int* nijk     = ipt_nijk     + 3*c_tg;
       E_Int* ijk_start= ipt_nijk     + 3*c_tg +3*mxzone;
       E_Int* ndimdx   = ipt_ndimdx   +   c_tg;

       E_Int* nozone   = ipt_nozone   +   c_tg;
       E_Int* nosszone = ipt_nosszone +   c_tg;

       E_Float* hpc_cups = ipt_HPC_CUPS+   c_tg;

       //printf("zone %d %d %d %d %d \n", c, ndimdx_max, c_tg, nozone[0],  nosszone[0]);

       E_Int* nijk_new     = ipt_nijk_new     + 3*c;
       E_Int* ijk_startnew = ipt_nijk_new     + 3*c +3*mxzone;
       E_Int* ndimdx_new   = ipt_ndimdx_new   +   c;
       E_Int* nozone_new   = ipt_nozone_new   +   c;
       E_Int* nosszone_new = ipt_nosszone_new +   c;

       E_Float* hpc_cups_new = ipt_HPC_CUPS_new +   c;

       nijk_new[0]     = nijk[0];
       nijk_new[1]     = nijk[1];
       nijk_new[2]     = nijk[2];
       ijk_startnew[0] = ijk_start[0];
       ijk_startnew[1] = ijk_start[1];
       ijk_startnew[2] = ijk_start[2];
       ndimdx_new[0]   = ndimdx[0];
       nozone_new[0]   = nozone[0];
       nosszone_new[0] = nosszone[0];
       ndimdx[0]       = -1;

       hpc_cups_new[0] = hpc_cups[0];
       // estimation des poids
       cout += nijk[0]*nijk[1]*nijk[2]/hpc_cups[0]; //cout en seconde
     } // loop zone

  E_Float cupsmoy = float(ndimt)/cout;
  //E_Float cupsmoy = CupsMax;
  E_Float poids;
  E_Float cells_tg_float=0;

  E_Int cells_tg = 0;
  E_Int flui_tg  = 0;
  E_Int fluj_tg  = 0;
  E_Int fluk_tg  = 0;
  
  FldArrayI inddm_zone(Nbzones*6);   E_Int* ipt_inddm_zone    = inddm_zone.begin();

  for (E_Int c = 0; c < Nbzones; c++)
     {  
        E_Int* nijk      = ipt_nijk_new     + 3*c;
        E_Int* ijk_start = ipt_nijk_new     + 3*c + 3*mxzone;
        E_Int* ind_dm    = ipt_inddm_zone   + 6*c;

        E_Int* nozone_new   = ipt_nozone_new   +   c;
        //E_Int  No_zone      = nozone_new[0];

        E_Float* hpc_cups = ipt_HPC_CUPS_new +   c;

        ind_dm[0]= ijk_start[0];
        ind_dm[2]= ijk_start[1];
        ind_dm[4]= ijk_start[2];
        ind_dm[1]= ijk_start[0] + nijk[0]-1;
        ind_dm[3]= ijk_start[1] + nijk[1]-1;
        ind_dm[5]= ijk_start[2] + nijk[2]-1;


	E_Int size_c = nijk[0]*nijk[1]*nijk[2];
	E_Int size_i =(nijk[0]+1)*nijk[1]*nijk[2];
	E_Int size_j =(nijk[1]+1)*nijk[0]*nijk[2];
	E_Int size_k = 0;
        if(nijk[2] != 1) size_k = (nijk[2]+1)*nijk[0]*nijk[1];

        //poids = cupsmoy/ipt_HPC_CUPS[c];
        //poids = cupsmoy/ipt_HPC_CUPS[No_zone];
        poids = cupsmoy/hpc_cups[0];
        //if(nstep==1) printf("poids %f %f %f %d %d \n", poids, cupsmoy,ipt_HPC_CUPS[No_zone], c, No_zone);

        cells_tg_float += size_c*poids;
        flui_tg        += size_i*poids;
        fluj_tg        += size_j*poids;
        fluk_tg        += size_k*poids;
     }  

    cells_tg  = int(cells_tg_float) / __NUMTHREADS__;
    flui_tg  /= __NUMTHREADS__;
    fluj_tg  /= __NUMTHREADS__;
    fluk_tg  /= __NUMTHREADS__;

   E_Float cc = 1.; E_Float ci = 0.; E_Float cj=0.; E_Float ck= 0.;
   //E_Float cc = 0.19; E_Float ci = 0.25; E_Float cj=0.28; E_Float ck= 0.28;
   E_Int cost = cells_tg*cc + flui_tg*ci + fluj_tg*cj +fluk_tg*ck;
   E_Int cells_tg_save = cost;

   // creation et initialisation remainder
   FldArrayI tab_remaind_sav(__NUMTHREADS__); E_Int* remaind_sav  = tab_remaind_sav.begin();
   FldArrayI tab_remaind_pos(__NUMTHREADS__); E_Int* remaind_pos  = tab_remaind_pos.begin();
   FldArrayI tab_remaind(__NUMTHREADS__);     E_Int* remaind  = tab_remaind.begin();
   FldArrayI tab_cells(__NUMTHREADS__);       E_Int* cells    = tab_cells.begin();
   FldArrayI tab_flu_i(__NUMTHREADS__);       E_Int* flu_i    = tab_flu_i.begin();
   FldArrayI tab_flu_j(__NUMTHREADS__);       E_Int* flu_j    = tab_flu_j.begin();
   FldArrayI tab_flu_k(__NUMTHREADS__);       E_Int* flu_k    = tab_flu_k.begin();
   FldArrayI tab_grain(__NUMTHREADS__);       E_Int* grain    = tab_grain.begin();

   FldArrayI tab_topo_lu(3);   E_Int* topo_lu  = tab_topo_lu.begin();
   FldArrayI tab_ind_dm_th(6); E_Int* ind_dm_th= tab_ind_dm_th.begin();
 
   FldArrayI tab_dim_i(__NUMTHREADS__);   E_Int* dim_i  = tab_dim_i.begin();
   FldArrayI tab_dim_j(__NUMTHREADS__);   E_Int* dim_j  = tab_dim_j.begin();
   FldArrayI tab_dim_k(__NUMTHREADS__);   E_Int* dim_k  = tab_dim_k.begin();

   FldArrayI tab_numa_socket(NBR_SOCKET);   E_Int* numa_socket  = tab_numa_socket.begin();

   E_Int dependence[__NUMTHREADS__];

   for (E_Int th = 0; th < __NUMTHREADS__ ; th++){ remaind[th] = cost; cells[th]=0; flu_i[th]=0;  flu_j[th]=0;  flu_k[th]=0; grain[th]=0; dependence[th]=-1;}


   E_Int PtiterTask = ipt_omp[nssiter + nstep-1];

   for (E_Int c = 0; c < Nbzones; c++)
     {  
        E_Int* ind_dm       = ipt_inddm_zone   + 6*c;
        E_Int* nijk         = ipt_nijk_new     + 3*c;
        E_Int* nozone_new   = ipt_nozone_new   +   c;
        E_Int* nosszone_new = ipt_nosszone_new +   c;

        E_Int  No_zone      = nozone_new[0];
        E_Int  No_sszone    = nosszone_new[0];
        //E_Int* ndimdx       = ipt_ndimdx_new   +   c;

        E_Float* hpc_cups = ipt_HPC_CUPS_new +   c;

        E_Int* ipt_nidom_loc = ipt_ind_dm[No_zone] + param_int[No_zone][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
        E_Int nb_subzone     = ipt_nidom_loc [nstep-1];   

        //printf("c Nozone %d %d %d %d %d %d %d %d %d \n", c, No_zone, ind_dm[0], ind_dm[1], ind_dm[2], ind_dm[3], ind_dm[4], ind_dm[5], ndimdx[0] );

        E_Int PtTask = PtiterTask +c*(6+ __NUMTHREADS__*7);

        E_Int  PtTaskomp    = PtTask + 6+ __NUMTHREADS__*7 ;

        //printf("Ptomp= %d , ptiter= %d , pt_subzone= %d , nstep= %d \n ",Ptomp, PtrIterOmp, PtZoneomp ,   nstep);
        E_Int check = PtTaskomp - shift_omp;
        if(check > mx_omp_size_int && boom ==0) {boom =1; check_size =check;}

        ipt_omp[PtTask  ] = No_zone;
        ipt_omp[PtTask+1] = No_sszone;

        //initilisation compteur Numa
        if(nstep==1) { for (E_Int socket = 0; socket < NBR_SOCKET; socket++) { numa_socket[socket] = 0;} }
 
        poids = cupsmoy/hpc_cups[0];
        //poids = cupsmoy/ipt_HPC_CUPS[c];
        //poids = cupsmoy/ipt_HPC_CUPS[No_zone];
        //printf("Nozone %d %d %f %d \n",c, No_zone, hpc_cups[0], nstep); 

        E_Int size_c = nijk[0]*nijk[1]*nijk[2]*poids;
        E_Int size_i =(nijk[0]+1)*nijk[1]*nijk[2]*poids;
        E_Int size_j =(nijk[1]+1)*nijk[0]*nijk[2]*poids;
        E_Int size_k = 0;
        if(nijk[2] != 1) size_k = (nijk[2]+1)*nijk[0]*nijk[1]*poids;

        E_Int cells_tmp = size_c;

        size_c = (size_c*cc + size_i*ci + size_j*cj + size_k*ck);

        E_Int Nthreads;
        //printf("so %d %d  %f %d \n ", size_c, nijk[0]*nijk[1]*nijk[2], poids, c);

        //
        //on determine les threads ou il reste le plus de place
        for (E_Int th1 = 0; th1 < __NUMTHREADS__ ; th1++) { remaind_sav[th1] =  remaind[th1];}
        for (E_Int th1 = 0; th1 < __NUMTHREADS__ ; th1++)
        { 
          E_Int rmax=-1e9;
          for (E_Int th = 0; th < __NUMTHREADS__ ; th++) { if(remaind[th] > rmax){ rmax = remaind[th]; remaind_pos[th1]=th;} }
          remaind[ remaind_pos[th1] ]=-1e9;
        }
        // on reinitialise les restes apres le tri par ordre decroissant
        for (E_Int th1 = 0; th1 < __NUMTHREADS__ ; th1++) { remaind[th1] =  remaind_sav[th1];}

        E_Int rmax =  remaind[ remaind_pos[0] ];

        E_Float      marge  = 1.005;
        if      (rmax<= cells_tg_save/8. ){marge  = 1.16;} 
        else if (rmax<= cells_tg_save/4  ){marge  = 1.06;} 
        else if (rmax<= cells_tg_save/2  ){marge  = 1.04;} 
        else if (rmax<= cells_tg_save/1.5){marge  = 1.035;} 

        //printf("size_c/poid= %d %f , remain max= %d , th_max= %d,  marge= %f, zone= %d, target= %d \n ", size_c, poids, rmax, remaind_pos[0], marge, No_zone, cells_tg_save);

        if(size_c <= rmax*marge || size_c<= 0.01*cells_tg_save)  // si place sur 1 thread ou zone tres petite, on reserve 1 thread pour eviter synchro sur micozone
        {
            //E_Int socket_tg = c%nb_socket;
	          E_Int th_current =  __NUMTHREADS__ -1 -c%__NUMTHREADS__;
            E_Int sens =-1;
            E_Int dirr = c/ __NUMTHREADS__;
            if (dirr%2==0) {sens=1; th_current=c%__NUMTHREADS__;}
            //printf("verif=  dirr= %d , sens= %d,  thcurrent= %d, zone= %d,  \n ", dirr, sens, th_current,  c);
            
	    E_Int th_count   = 0;
	    //E_Int th_current = socket_tg*core_per_socket;

            E_Int lgo = 1;
	    while(lgo==1)
              { //printf("th_current mono th %d %d \n", lgo, c);

                if(size_c/marge <=  remaind[th_current]) lgo = 0;

                //on ajuste le sens de parcours pour la recherche de place si on est au limite min ou max des threads
                if (th_current==0                && sens==-1) sens= 1;
                if (th_current==__NUMTHREADS__-1 && sens== 1) sens=-1;

                if(lgo==1 && th_count == __NUMTHREADS__-1)
                   {
                     lgo = 0;
                     //on cherche le thread le moins charge pour attribuer le residu
                     E_Int th_max =-100000;
                     for (E_Int th = 0; th < __NUMTHREADS__; th++)
                       {  
                        if(remaind[th] > th_max){ th_max = remaind[th]; th_current = th;}
                       }
                   }
                if(lgo==1)
                { th_current +=sens;
                  th_count   +=1;
                  if(th_current==__NUMTHREADS__-1 && sens==1) th_current = 0;
                  if(th_current== 0 && sens==-1) th_current = __NUMTHREADS__-1;
                }
              }// while go
  
            remaind[th_current] =  remaind[th_current] - size_c; 

            cells[th_current] +=  cells_tmp;
            flu_i[th_current] +=  size_i; 
            flu_j[th_current] +=  size_j; 
            flu_k[th_current] +=  size_k; 

            grain[th_current] +=1;
             
            //Carte threads actifs       
            for (E_Int th = 0; th < __NUMTHREADS__; th++) { ipt_omp[PtTask + 2 + th] = -2; } //Thread inactif

            //le Thread "0" (pas au sens omp: indirection entre omp_numthread et ithread)) prend le job
            ipt_omp[PtTask+ 2 + th_current] = 0;
            //Nb threads actifs 
            Nthreads = 1;
            ipt_omp[PtTask + __NUMTHREADS__ +2 ] = Nthreads ;
            //topology omp
            ipt_omp[PtTask + __NUMTHREADS__ +3 ] = 1 ;
            ipt_omp[PtTask + __NUMTHREADS__ +4 ] = 1 ;
            ipt_omp[PtTask + __NUMTHREADS__ +5 ] = 1 ;
            //indice sous-domaine omp

            ipt_omp[PtTask + __NUMTHREADS__ +6 ] =  ind_dm[0];
            ipt_omp[PtTask + __NUMTHREADS__ +7 ] =  ind_dm[1];
            ipt_omp[PtTask + __NUMTHREADS__ +8 ] =  ind_dm[2];
            ipt_omp[PtTask + __NUMTHREADS__ +9 ] =  ind_dm[3];
            ipt_omp[PtTask + __NUMTHREADS__ +10] =  ind_dm[4];
            ipt_omp[PtTask + __NUMTHREADS__ +11] =  ind_dm[5];

	   if ((display==1 && nitrun ==1 && nstep==1) || display ==2)
           {
             printf("souszone %d de la zone %d calcule par %d  Threads,  nstep= %d, cycle = %d \n", c, No_zone, Nthreads, nstep, param_int[No_zone][NSSITER]/param_int[No_zone][LEVEL] );
             printf("-- Thread= %d : fenetre=( %d, %d, %d, %d, %d, %d) \n", th_current, ind_dm[0],ind_dm[1],ind_dm[2],ind_dm[3],ind_dm[4],ind_dm[5]);
           }
        }
        else
        {
            Nthreads=0; E_Int cells_loc=0; E_Int th_min; E_Int cells_tg_min; E_Int cells_tg_loc;
            // zone plus grande que la cible globale
            if(size_c*3 > cells_tg_save)
            {
              while( cells_loc < size_c && Nthreads < __NUMTHREADS__ )
	        {  
	    	  //cells_loc += remaind[ remaind_pos[Nthreads]]; Nthreads +=1;
		      cells_loc += remaind[ remaind_pos[Nthreads]]*marge; Nthreads +=1; // modif avril

		      //printf("cells_loc = %d, remain = %d %d %d\n", cells_loc, remaind[ remaind_pos[Nthreads-1]], Nthreads, size_c ); 
	        }

                cells_tg_loc  =  cells_loc/Nthreads;


                cells_tg_min = remaind[ remaind_pos[Nthreads-1]]; th_min=-1; 
           
                if( cells_tg_min != cells_tg_save){ th_min=remaind_pos[Nthreads-1]; } 
                if( cells_tg_loc == cells_tg_min) th_min =-1;

                cells_tg_loc  =  cells_tg_min*marge; //modif avril
                cells_tg_loc  = K_FUNC::E_max( cells_tg_loc , rmax );
	    }
            else
            //{     E_Int target = cells_tg_save*(marge-1); 
            {     E_Int target = cells_tg_save*(1.06-1); 
                  Nthreads     = size_c/target+1;
                  if(Nthreads > __NUMTHREADS__) Nthreads = __NUMTHREADS__;
                  cells_loc    = size_c/Nthreads;
                  cells_tg_loc =  cells_loc;
                  cells_tg_min = cells_loc; th_min=-1; 
            }

            //printf("cell_tg %d, rmax= %d , cell_thmin= %d, th_min= %d , Nthreads= %d , zone %d \n", cells_tg, rmax, cells_tg_min, th_min,Nthreads, c );
            //if(c==8) printf("cell_tg %d %d , th_min= %d \n", cells_tg, rmax, th_min);
          E_Int nd =No_zone;
#include "FastC/HPC_LAYER/SIZE_MIN.cpp"


            E_Int search      = 1;
            E_Int search_topo = 0;
            E_Int adapt_thread= 0;
            while(search == 1)
              {
               //printf("th_current multi th %d %d %d \n", search, Nthreads, c);
               //
                //verif place dispo sur les dernier threads
                //for th_check in range(OMP_NUM_THREADS-1,OMP_NUM_THREADS-Nthreads-1,-1):
                //   print 'remain', remaind[str(th_check)], th_check

                search = 0;

                //printf("Nthread %d \n", Nthreads);
                E_Int ithread = 1;
                if(search_topo==0)
                      {
                        indice_boucle_lu_(c, ithread, Nthreads, lmin, ind_dm, topo_lu, ind_dm_th );
                        E_Int Nthread_tg = topo_lu[0]*topo_lu[1]*topo_lu[2];

                        while(Nthread_tg == 1 && Nthreads !=1)
                         { Nthreads-=1;
                          indice_boucle_lu_(c, ithread, Nthreads, lmin, ind_dm, topo_lu, ind_dm_th );
                          Nthread_tg = topo_lu[0]*topo_lu[1]*topo_lu[2];
                         } 
                      }
                //if(c==8) printf("nijk      a %d %d %d  \n", nijk[0], nijk[1],nijk[2]);
                //printf("topo_LLLUUU %d %d %d %d \n", topo_lu[0], topo_lu[1],topo_lu[2], search_topo);
                //if(c==0) printf("topo_LLLUUU %d %d %d %d \n", topo_lu[0], topo_lu[1],topo_lu[2], search_topo);

                
                for (E_Int dir = 0; dir < 3; dir++){
                  for (E_Int l = 0; l < topo_lu[dir]; l++){

                        if     (dir==0){dim_i[l] =  nijk[dir]/topo_lu[dir]; if(l < nijk[dir]%topo_lu[dir]){dim_i[l] += 1;}  }
                        else if(dir==1){dim_j[l] =  nijk[dir]/topo_lu[dir]; if(l < nijk[dir]%topo_lu[dir]){dim_j[l] += 1;}  }
                        else           {dim_k[l] =  nijk[dir]/topo_lu[dir]; if(l < nijk[dir]%topo_lu[dir]){dim_k[l] += 1;}  }
                    }
                  }

                /*for (E_Int dir = 0; dir < 3; dir++){
                  for (E_Int l = 0; l < topo_lu[dir]; l++){

                        if     (dir==0){printf("dim i %d  \n", dim_i[l]);  }
                        else if(dir==1){printf("dim j %d  \n", dim_j[l]);  }
                        else           {printf("dim k %d  \n", dim_k[l]);  }
                    }
                  }
                */

                //poids = cupsmoy/ipt_HPC_CUPS[No_zone];
                poids = cupsmoy/hpc_cups[0];
                //poids = cupsmoy/ipt_HPC_CUPS[c];

                E_Int res = (cc*dim_i[0]*dim_j[0]*dim_k[0] + ci*dim_i[0]*dim_j[0]*dim_k[0] + cj*dim_j[0]*dim_i[0]*dim_k[0] + ck*dim_k[0]*dim_i[0]*dim_j[0] )*poids - cells_tg_loc;
                E_Float sign = 1.;
                if(res < 0) sign = -1.; //block trop petit % la cible

                //#print 'size bloc=',dims_i[0]*dim_j[0]*dim_k[0 , 'sizetg=', cells_tg_loc
            
                //sous bloc trop grand
                if(res*sign > 0 && adapt_thread ==0)
                  {
                    //E_Int compteur =0;
                    E_Int go = 0;
                    //while(res*sign > 0. && go==0 && compteur < compteur_mx)
                    while(res*sign > 0. && go==0)
                       {
                        go      = 1;
                        E_Int deriv_i= K_FUNC::E_min(1, (topo_lu[0]-1) );
                        E_Int deriv_j= K_FUNC::E_min(1, (topo_lu[1]-1) );
                        E_Int deriv_k= K_FUNC::E_min(1, (topo_lu[2]-1) );

                        E_Int cout_i = dim_j[0]*dim_k[0]*deriv_i;
                        E_Int cout_j = dim_i[0]*dim_k[0]*deriv_j;
                        E_Int cout_k = dim_j[0]*dim_i[0]*deriv_k;
                        if( (cout_i+cout_j+cout_k)*poids < res*sign)
                          {
                            dim_i[0] -=deriv_i*sign; 
                            dim_j[0] -=deriv_j*sign;
                            dim_k[0] -=deriv_k*sign;
                            res      -= (cout_i+cout_j+cout_k)*sign*poids;
                            if(cout_i+cout_j+cout_k != 0) go = 0;
                //if(c==8) printf("cout  %d %d %d  \n", res, cout_k,dim_k[0]);
                          }
                        else if( (cout_j+cout_k)*poids < res*sign)
                          {
                            dim_j[0] -=deriv_j*sign;
                            dim_k[0] -=deriv_k*sign;
                                 res -= (cout_j+cout_k)*sign*poids;
                            if(cout_j+cout_k != 0) go = 0;
                          }
                        else if( (cout_j+cout_i)*poids < res*sign)
                          {
                            dim_j[0] -=deriv_j*sign; 
                            dim_i[0] -=deriv_i*sign;
                            res      -= (cout_j+cout_i)*sign*poids;
                            if(cout_j+cout_i != 0) go = 0;
                          }
                        else if( (cout_k+cout_i)*poids < res*sign)
                          {
                            dim_k[0] -=deriv_k*sign; 
                            dim_i[0] -=deriv_i*sign;
                            res      -= (cout_k+cout_i)*sign*poids;
                            if( cout_k+cout_i != 0) go = 0;
                          }
                        else if(cout_k*poids < res*sign)
                          {
                            dim_k[0] -=deriv_k*sign;
                            res      -= cout_k*sign*poids;
                            if(cout_k != 0) go = 0;
                          }
                        else if(cout_j*poids < res*sign)
                          {
                            dim_j[0] -=deriv_j*sign;
                            res      -= cout_j*sign*poids;
                            if(cout_j != 0) go = 0;
                          }
                        else if(cout_i*poids < res*sign)
                          {
                            dim_i[0] -=deriv_i*sign;
                            res      -= cout_i*sign*poids;
                            if(cout_i != 0) go = 0;
                          }
                        E_Int rest = nijk[0] - deriv_i*dim_i[0];
                        if(rest <= lmin || dim_i[0] <= lmin)  go=1;
                        rest = nijk[1] - deriv_j*dim_j[0];
                        if(rest <= lmin || dim_j[0] <= lmin)  go=1;
                        rest = nijk[2] - deriv_k*dim_k[0];
                        //if(c==0) printf("rest %d %d  \n",dim_k[0], lmin );
                        if(rest <= lmin || dim_k[0] <= lmin)  go=1;
                        //compteur +=1;
                       }//while go
                  } // if res

              //printf("topo dim0 %d %d %d  \n",dim_i[0], dim_j[0],dim_k[0] );
              //if(c==0) printf("topo dim0 %d %d %d  \n",dim_i[0], dim_j[0],dim_k[0] );

              for (E_Int dir = 0; dir < 3; dir++){
                   if(topo_lu[dir] != 1)
                   { 
                      //on repartie la place restante entre les thread restant
                      if(th_min != -1)
                      {

                          E_Int sum;
                          if(topo_lu[0]*topo_lu[1]==1)  //decoupe K
                           {
                            sum=0;
                            for (E_Int l = 0; l <topo_lu[2] ; l++)
                            {
		              E_Float flagk= 0; if(l== topo_lu[2]-1) flagk+=ck;
                              dim_k[l]=(remaind[ remaind_pos[l]] -flagk*nijk[0]*nijk[1]*poids) /( nijk[0]*nijk[1]*cc + ci*(nijk[0]+1)*nijk[1] + cj*(nijk[1]+1)*nijk[0] + ck*nijk[1]*nijk[0] )/poids;
                              dim_k[l]= K_FUNC::E_max(dim_k[l], lmin);
                              if (l==topo_lu[2]-1){dim_k[l] = nijk[2]-sum;}
                              sum+=dim_k[l];
                            }
                           }
                          else if (topo_lu[0]*topo_lu[2]==1) //decoupe J
                           {
                            sum=0;
                            for (E_Int l = 0; l <topo_lu[1] ; l++)
                            {
		              E_Float flagj= 0; if(l== topo_lu[1]-1) flagj+=cj;
                              dim_j[l]=(remaind[ remaind_pos[l]] -flagj*nijk[0]*nijk[2]*poids) /( nijk[0]*nijk[2]*cc + ci*(nijk[0]+1)*nijk[2] + ck*(nijk[2]+1)*nijk[0] + cj*nijk[2]*nijk[0] )/poids;
                              dim_j[l]= K_FUNC::E_max(dim_j[l], lmin);
                              if (l==topo_lu[1]-1){dim_j[l] = nijk[1]-sum;}
                              sum+=dim_j[l];
                            }
                           }
                          else if (topo_lu[1]*topo_lu[2]==1) //decoupe I
                           {
                            sum=0;
                            for (E_Int l = 0; l <topo_lu[0] ; l++)
                            {
		              E_Float flagi= 0; if(l== topo_lu[0]-1) flagi+=ci;
                              dim_i[l]=(remaind[ remaind_pos[l]] -flagi*nijk[1]*nijk[2]*poids) /( nijk[1]*nijk[2]*cc + cj*(nijk[1]+1)*nijk[2] + ck*(nijk[2]+1)*nijk[1] + ci*nijk[2]*nijk[1] )/poids;
                              dim_i[l]= K_FUNC::E_max(dim_i[l], lmin);
                              if (l==topo_lu[0]-1){dim_i[l] = nijk[0]-sum;}
                              sum+=dim_i[l];
                            }
                           }
                          else if (topo_lu[2]==1) //decoupe IJ
                           {

			     //printf("DECOUPE IJ \n");

                            sum=0; // on dimensione dim_i a partie de dim_j(0)
                            for (E_Int l = 0; l <topo_lu[0] ; l++)
                            {
                              dim_i[l]= remaind[ remaind_pos[l] ]/ ( dim_j[0]*nijk[2]*( cc + cj + ci) + ck*(nijk[2]+1)*dim_j[0]  )/poids;
                              //
                              dim_i[l]= K_FUNC::E_max(dim_i[l], lmin);
                              if (l==topo_lu[0]-1){dim_i[l] = nijk[0]-sum;}
                              //printf("dimI %d %d \n", dim_i[l], l);
                              sum+=dim_i[l];
                            }
                            E_Int sum=dim_j[0]; // on dimensione dim_j a partie de dim_i(0)

                            for (E_Int l = 1; l <topo_lu[1] ; l++)
                            {
                              dim_j[l]=remaind[ remaind_pos[l*topo_lu[0]]] /( dim_i[0]*nijk[2]*(cc + ci +cj) + ck*(nijk[2]+1)*dim_i[0]  )/poids;
                              //
                              dim_j[l]= K_FUNC::E_max(dim_j[l], lmin);
                              if (l==topo_lu[1]-1){dim_j[l] = nijk[1]-sum;}
                              sum+=dim_j[l];
                              //printf("dimJ %d %d \n", dim_j[l],dim_j[0]);
                            }
                           }
                          else if (topo_lu[0]==1) //decoupe KJ
                           {
                            sum=0; // on dimensione dim_k a partie de dim_j(0)


                            for (E_Int l = 0; l <topo_lu[2] ; l++)
                            {
                              dim_k[l]= remaind[ remaind_pos[l]] /( nijk[0]*dim_j[0]*( cc +cj +ck) + ci*(nijk[0]+1)*dim_j[0] )/poids;
                              dim_k[l]= K_FUNC::E_max(dim_k[l], lmin);
                              if (l==topo_lu[2]-1){dim_k[l] = nijk[2]-sum;}
                              //printf("dimK %d %d \n", dim_k[l], l);
                              sum+=dim_k[l];
                            }
                            E_Int sum=dim_j[0]; // on dimensione dim_j a partie de dim_k(0)

                            for (E_Int l = 1; l <topo_lu[1] ; l++)
                            {
                              dim_j[l]= remaind[ remaind_pos[ l*topo_lu[2]]] /( dim_k[0]*nijk[0]*(cc + cj +ck) + ci*(nijk[0]+1)*dim_k[0] )/poids;
                              dim_j[l]= K_FUNC::E_max(dim_j[l], lmin);
                              if (l==topo_lu[1]-1){dim_j[l] = nijk[1]-sum;}
                              //printf("dimJ %d %d \n", dim_j[l], l);
                              sum+=dim_j[l];
                            }
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

              //En priorite, on cherche si un changement de topo est possible
              if(search_topo == 0)
              {
                 E_Int change_topo=0;
                 for (E_Int dir = 0; dir < 3; dir++){
                  for( E_Int l=0; l < topo_lu[dir]; l++){
                       E_Int test;
                       if     (dir==0) test = dim_i[l];
                       else if(dir==1) test = dim_j[l];
                       else            test = dim_k[l];
                       if(search==0 && test < lmin &&  nijk[dir] != test) {change_topo=1;}  
                     }//loop topo
                   } // loop dir

                 if (change_topo==1)
                 { 
                  search_topo = topo_test( topo_lu, nijk, cells_tg_loc, lmin, dim_i[0], dim_j[0], dim_k[0] );
                  search = 1;
                  //printf("change topo %d %d %d %d \n", topo_lu[0], topo_lu[1],topo_lu[2], search_topo);
                 }
                 //search_topo = 1;  
              }
              //Sinon on reduit le nb de threads
              else
              {
               for (E_Int dir = 0; dir < 3; dir++){
                for( E_Int l=0; l < topo_lu[dir]; l++){
                      E_Int test;
                      if     (dir==0) test = dim_i[l];
                      else if(dir==1) test = dim_j[l];
                      else            test = dim_k[l];
                      if(search==0 && test < lmin &&  nijk[dir] != test) 
                      {
                          Nthreads -=  1;

                          search = 1;
                          adapt_thread =1;
                          //printf("change th   %d \n", Nthreads);
                      }
                    }//loop topo
                  } // loop dir
               search_topo = 0;  
               }
              }// while search



            E_Int list_affected[__NUMTHREADS__];
            E_Int   thread_list[__NUMTHREADS__];
            //Carte threads actifs       
            for (E_Int th = 0; th < __NUMTHREADS__; th++)
              { list_affected[th]=-1; 
                ipt_omp[PtTask + 2 + th] = -2;} //Thread inactif

            //poids = cupsmoy/ipt_HPC_CUPS[No_zone];
            poids = cupsmoy/hpc_cups[0];
            //poids = cupsmoy/ipt_HPC_CUPS[c];
            //Recherche thread disponible
            E_Int th = 0;
            for (E_Int k = 0; k < topo_lu[2]; k++){
              for (E_Int j = 0; j < topo_lu[1]; j++){
                for (E_Int i = 0; i < topo_lu[0]; i++){

                  //E_Int socket_tg  = c%nb_socket;
                  E_Int th_count   = 0;
                  E_Int th_current =  __NUMTHREADS__ -1 -c%__NUMTHREADS__;
                  E_Int sens =-1;
                  E_Int dirr = c/ __NUMTHREADS__;
                  if (dirr%2==0) {sens=1; th_current=c%__NUMTHREADS__;}
	          //E_Int th_current = socket_tg*core_per_socket;

		  //E_Int th_current = 0;
		  E_Int szi= dim_i[i];
		  E_Int szj= dim_j[j];
		  E_Int szk= dim_k[k];
                  if(i== topo_lu[0]-1) szi+=1;
                  if(j== topo_lu[1]-1) szj+=1;
                  if(k== topo_lu[2]-1) szk+=1;
                  E_Int size_th = (cc*dim_i[i]*dim_j[j]*dim_k[k] + ci*szi*dim_j[j]*dim_k[k] + cj*szj*dim_i[i]*dim_k[k] + ck*szk*dim_i[i]*dim_j[j] )*poids;
                  E_Int size_loc   = size_th;

	          //On cherche de la place sur un des socket, puis sur le suivant si pas de candidat local
	          while( (size_loc/marge > remaind[th_current] || list_affected[th_current]!=-1)  && th_count < __NUMTHREADS__-1) 
                    { 
                      //on ajuste le sens de parcours pour la recherche de place si on est au limite min ou max des threads
                      if (th_current==0                && sens==-1) sens= 1;
                      if (th_current==__NUMTHREADS__-1 && sens== 1) sens=-1;

                      th_current +=sens; 
                      th_count   +=1;
                      if(th_current==__NUMTHREADS__-1 && sens==1) th_current = 0;
                      if(th_current== 0 && sens==-1) th_current = __NUMTHREADS__-1;
                      //if(th_current==__NUMTHREADS__-1) th_current = 0;
                    }
                
                  if(th_count == __NUMTHREADS__-1 || list_affected[th_current]!=-1)
                  {
                     E_Int blind = -1000000000;
                     for (E_Int l = 0; l < __NUMTHREADS__; l++){
                        if(remaind[l] > blind && list_affected[l]==-1)
                        { blind = remaind[l];
                          th_current = l;
                        }
                     }
                  }
                  //tentative optim distrib
                  //if(th_min != -1) { if(list_affected[th_min]==-1){ th_current = th_min; th_min = -1;} }   

                  //on exclut le thread our ne pas avoir 2 souszone traite par le meme thread
                  list_affected[th_current]=0;

                  //on estime a quelle socket Numa la zone appartient  ( a blinder...)
                  E_Int socket = th_current/core_per_socket;

                  if (socket >= nb_socket || socket >= NBR_SOCKET) socket=0;
                  numa_socket[socket] += size_loc;

                  remaind[th_current] -=  size_loc;
                  E_Int size_i =  dim_i[i]*dim_j[j]*dim_k[k];
                  E_Int size_j =  dim_i[i]*dim_j[j]*dim_k[k];
                  E_Int size_k =  dim_i[i]*dim_j[j]*dim_k[k];
                  if(i == topo_lu[0]-1) { size_i = (dim_i[i]+1)*dim_j[j]*dim_k[k];}
                  if(j == topo_lu[1]-1) { size_j = dim_i[i]*(dim_j[j]+1)*dim_k[k];}
                  if(k == topo_lu[2]-1) { size_k = dim_i[i]*dim_j[j]*(dim_k[k]+1);}
                  if(nijk[2] == 1) size_k = 0;
                
                  cells[th_current]   +=  dim_i[i]*dim_j[j]*dim_k[k];
                  flu_i[th_current]   +=  size_i*poids; 
                  flu_j[th_current]   +=  size_j*poids; 
                  flu_k[th_current]   +=  size_k*poids; 

                  grain[th_current] +=1;

                  //Thread actif
                  ipt_omp[PtTask + 2 + th_current] = th;  
                  thread_list[th] = th_current;

		  //printf("No_zone,  PtZoneomp + th_current, th = %d %d %d %d\n",No_zone, PtZoneomp + th_current, grain[th_current], th_current); 

                  //pour optimiser graph dependance, on bascule en premiere position les threads deja sollicité sur les zones precedente
		  
                  if( dependence[th_current]==1)
                  {  
                    E_Int thcurrent_save, size;
                    if ( i == topo_lu[0]-1 && i !=0)
                     { size    = dim_i[i];
                       dim_i[i]   = dim_i[0]; 
                       dim_i[0]   = size;
                       thcurrent_save = thread_list[th-i];
                       th_current = thread_list[th];
                       ipt_omp[PtTask + 2 + thcurrent_save] = th;  
                       ipt_omp[PtTask + 2 + th_current] = th-i;  
		       thread_list[th-i] = th_current;
		       thread_list[th] = thcurrent_save;		       
                     }
                    if ( j == topo_lu[1]-1  && j !=0)
                     { size    = dim_j[j];
                       dim_j[j]= dim_j[0]; 
                       dim_j[0]= size;
                       thcurrent_save = thread_list[th-j];
                       th_current = thread_list[th];
                       ipt_omp[PtTask + 2 + thcurrent_save] = th;  
                       ipt_omp[PtTask + 2 + th_current] = th-j;  
		       thread_list[th-j] = th_current;
		       thread_list[th] = thcurrent_save;
		       //printf("th, th-j, thcurrent_save, th_current = %d %d %d %d \n", th, th-j, thcurrent_save, th_current);
                     }
                    if ( k == topo_lu[2]-1  && k !=0)
                     { size    = dim_k[k];
                       dim_k[k]= dim_k[0]; 
                       dim_k[0]= size;
                       thcurrent_save = thread_list[th-k];
                       th_current = thread_list[th];
                       ipt_omp[PtTask + 2 + thcurrent_save] = th;  
                       ipt_omp[PtTask + 2 + th_current] = th-k;  
		       thread_list[th-k] = th_current;
		       thread_list[th] = thcurrent_save;
		       //printf("th, th-k, thcurrent_save, th_current = %d %d %d %d \n", th, th-k, thcurrent_save, th_current);
                     }
                  }
		  
                  if(i == topo_lu[0]-1) { dependence[th_current]=1;}
                  if(j == topo_lu[1]-1) { dependence[th_current]=1;}
                  if(k == topo_lu[2]-1) { dependence[th_current]=1;}

                  //printf("th_currentNew %d %d %d %d %d %d %d %d \n", th_current, th, remaind[th_current], i,j,k, size_loc, nstep  );

                  th +=1;

                  //print 'thread=', th, size_loc, remaind[str(th_current)], th_current
                }//loop i
              }//loop j
            }//loop k

            if(nstep ==1)
            {E_Int numa_max=0; E_Int sock_tg =0;
             for (E_Int socket = 0; socket < NBR_SOCKET; socket++) { if( numa_socket[socket] > numa_max) { numa_max= numa_socket[socket]; sock_tg= socket;} } 
             //Ivan param_int[No_zone][Ptomp + nssiter]= sock_tg;
            }

            
            //Nb threads actifs
            ipt_omp[PtTask + __NUMTHREADS__ +2 ] =  topo_lu[0]*topo_lu[1]*topo_lu[2];
            //topology omp
            ipt_omp[PtTask + __NUMTHREADS__ +3 ] =  topo_lu[0];
            ipt_omp[PtTask + __NUMTHREADS__ +4 ] =  topo_lu[1];
            ipt_omp[PtTask + __NUMTHREADS__ +5 ] =  topo_lu[2];

            th = 0;
            E_Int kstart = ind_dm[4];
            for (E_Int k = 0; k < topo_lu[2]; k++){
              ind_dm_th[4] = kstart;
              ind_dm_th[5] = kstart + dim_k[k] - 1;

              E_Int jstart = ind_dm[2];
              for (E_Int j = 0; j < topo_lu[1]; j++){
                ind_dm_th[2] = jstart;
                ind_dm_th[3] = jstart + dim_j[j] - 1;

                E_Int istart = ind_dm[0];
                for (E_Int i = 0; i < topo_lu[0]; i++){


                    ind_dm_th[0] = istart;
                    ind_dm_th[1] = istart + dim_i[i] - 1;
                    istart       = ind_dm_th[1] + 1;
                    
                    ipt_omp[PtTask + __NUMTHREADS__ +6 + th*6] =  ind_dm_th[0];
                    ipt_omp[PtTask + __NUMTHREADS__ +7 + th*6] =  ind_dm_th[1];
                    ipt_omp[PtTask + __NUMTHREADS__ +8 + th*6] =  ind_dm_th[2];
                    ipt_omp[PtTask + __NUMTHREADS__ +9 + th*6] =  ind_dm_th[3];
                    ipt_omp[PtTask + __NUMTHREADS__ +10+ th*6] =  ind_dm_th[4];
                    ipt_omp[PtTask + __NUMTHREADS__ +11+ th*6] =  ind_dm_th[5];
		    //printf(" No_zone, adresse = %d  %d \n", No_zone,  PtZoneomp + __NUMTHREADS__ +4 + th*6);
                    //printf(" dom = %d  %d  %d  %d %d %d %d \n", ind_dm_th[0],ind_dm_th[1], ind_dm_th[2],ind_dm_th[3],ind_dm_th[4], ind_dm_th[5], th);
                    th     +=1;
                }//loop i
               jstart  =   ind_dm_th[3]+ 1;
              }//loop j
             kstart  = ind_dm_th[5] + 1;
            }//loop k

            E_Int  PtTaskomp  = PtTask + 6+ __NUMTHREADS__*7 ;
            E_Int check       = PtTaskomp - shift_omp;
            if(check > mx_omp_size_int && boom ==0) {boom =1; check_size =check;}

           if ((display==1 && nitrun ==1 ) || display==2)
           { 
             E_Int Nthreads = ipt_omp[PtTask + __NUMTHREADS__ +2 ];

             printf("souszone %d de la zone %d calcule par %d  Threads,  nstep= %d, cycle=%d \n", c, No_zone, Nthreads, nstep, param_int[No_zone][NSSITER]/param_int[No_zone][LEVEL]);
             //printf("souszone %d de la zone %d calcule par %d  Threads,  nstep= %d, cycle=%d %d  \n", c, No_zone, Nthreads, nstep, PtTask,PtiterTask );
	     
             for (E_Int th_current = 0; th_current < __NUMTHREADS__; th_current++){

                 E_Int th = ipt_omp[PtTask + 2 + th_current] ; //Thread inactif
		 //printf("No_zone,  PtZoneomp + th_current, th = %d %d %d \n",No_zone, PtZoneomp + th_current, th); 

		 if( th != -2)
                 {
                   E_Int* ind_dm = ipt_omp + PtTask + __NUMTHREADS__ +6 + th*6;

                   printf("-- Thread= %d  %d : fenetre=( %d, %d, %d, %d, %d, %d) \n", th_current, th,
                                                 ind_dm[0],ind_dm[1],ind_dm[2],ind_dm[3],ind_dm[4],ind_dm[5]);
                 }
             }//loop
           }//if display

        }// if mono ou subzone

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
                printf("Fast.FastC.MX_OMP_SIZE_INT= %d\n ", check_size+20);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }

     }//loop zone

   //mise a jour Pointer pour stockage nstep+1
   if(nstep != nssiter)
   {
     ipt_omp[ nssiter +nstep  ]= ipt_omp[ nssiter +nstep-1  ] + Nbzones*(6+7*__NUMTHREADS__) ;
     //printf("ptiter distriv %d %d %d %d %d \n", ipt_omp[ nssiter +nstep],  ipt_omp[ nssiter +nstep-1], Nbzones, nstep-1, nitrun );

     for (E_Int c = 0; c < Nbzones; c++)
     {  
        E_Int* nozone_new   = ipt_nozone_new   +   c;
        E_Int  No_zone      = nozone_new[0];

        E_Int* ipt_nidom_loc = ipt_ind_dm[No_zone] + param_int[No_zone][ MXSSDOM_LU ]*6*nssiter + nssiter;     //nidom_loc(nssiter)
        E_Int nb_subzone     = ipt_nidom_loc [nstep-1];


        E_Int  PtTaskomp  = ipt_omp[ nssiter +nstep  ];
        E_Int check_size  = PtTaskomp - shift_omp;
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
                printf("Fast.FastC.MX_OMP_SIZE_INT= %d\n ", check_size+20);
                printf("------\n");
                printf("End error msg\n");
                printf("------\n");
                exit(0);
               }
     }  
   }

  if (display==1 || display ==2)
  //if ( (display==1 && nitrun ==1 && nstep == 1) || display ==2)
  //if ( (display==1 && nitrun ==1 ) || display ==2)
  {
     printf("OMP balancing: positive values --> slow thread. Target cells per Thread: %d Nitrun= %d, nstep= %d, cycle=%d \n", cells_tg_save, nitrun , nstep, lev);
     E_Float sum_c, sum_i,sum_j,sum_k;
     for (E_Int th_current = 0; th_current < __NUMTHREADS__; th_current++){
        if(fluk_tg!=0)
        {   E_Float dc =  float(cells[th_current]-cells_tg)/cells_tg*100;
            //E_Float dc =  (float(-remaind[th_current])/cells_tg_save)*100.;

            E_Float di =  float(flu_i[th_current]-flui_tg)/flui_tg*100;
            E_Float dj =  float(flu_j[th_current]-fluj_tg)/fluj_tg*100;
            E_Float dk =  float(flu_k[th_current]-fluk_tg)/fluk_tg*100;
            sum_c+= dc;sum_i+= di;sum_j+= dj;sum_k+= dk;
            printf("desequilibre thread, %d, : %f percent. Nb flui= %f,   Nb fluj= %f,  Nb fluk= %f, Nb zones= %d  \n", th_current, (float(-remaind[th_current])/cells_tg_save)*100., di, dj, dk,grain[th_current] );
        }
        else
        {  printf("desequilibre thread, %d, : %f percent. Nb flui= %f,   Nb fluj= %f,  Nb zones= %d \n", th_current, (float(-remaind[th_current])/cells_tg_save)*100., 
            float(flu_i[th_current]-flui_tg)/flui_tg*100, float(flu_j[th_current]-fluj_tg)/fluj_tg*100, grain[th_current] );
        }
     }
    //printf("verif globale, : %f  percent. Nb flui= %f,   Nb fluj= %f,  Nb fluk= %f \n", sum_c, sum_i,sum_j,sum_k );
  }
////
//
//
 } //fin loop sur les niveau en temps
///
//
//
  return;
}


// -----------------------------------------------------------------------------------
// 
// distributeThreads (python layer)
// 
// -----------------------------------------------------------------------------------
PyObject* K_FASTC::distributeThreads(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* dtloc;
  E_Int nssiter; E_Int nstep; E_Int nitrun; E_Int mx_omp_size_int; E_Int display;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOOlllll", &zones , &metrics, &dtloc, &nstep,  &nitrun, &mx_omp_size_int,  &display)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "OOOiiiii", &zones , &metrics, &dtloc, &nstep,  &nitrun, &mx_omp_size_int,  &display)) return NULL; 
#endif
  
  E_Int nidom = PyList_Size(zones);

  E_Int**   ipt_param_int;  E_Int** ipt_ind_dm; E_Float** ipt_param_real;
  ipt_param_int     = new E_Int*[nidom*2];
  ipt_ind_dm        = ipt_param_int   + nidom;
  ipt_param_real    = new E_Float*[nidom];


  vector<PyArrayObject*> hook;

  E_Int* ipt_dtloc =  K_PYTREE::getValueAI(dtloc, hook);

  for (E_Int nd = 0; nd < nidom; nd++)
  {    
    PyObject* zone = PyList_GetItem(zones, nd);
       
    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);
                       o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(o, hook);

    // get metric
    PyObject* metric     = PyList_GetItem(metrics, nd); // metric du domaine i
    ipt_ind_dm[nd]       = K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );
  }

  distributeThreads_c( ipt_param_int ,  ipt_param_real, ipt_ind_dm, nidom, ipt_dtloc, mx_omp_size_int, nstep, nitrun, display );

  delete [] ipt_param_int;
  delete [] ipt_param_real;

  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;
}

// -----------------------------------------------------------------------------------
// 
// optimisation topo thraed pour balancing
// 
// -----------------------------------------------------------------------------------
E_Int K_FASTC::topo_test( E_Int* topo, E_Int* nijk, E_Int& cells_tg, E_Int& lmin, E_Int& dim_i,  E_Int& dim_j, E_Int& dim_k)
{
  E_Int test =1;
   E_Int dim_loc, res;
  //
  // decoupe K
  //
  if (topo[0]*topo[1] == 1) 
  { //on test si on peut intervertir la direction K vers J
    dim_loc = cells_tg/(nijk[0]*nijk[2]);
    res     = nijk[1]-(topo[2]-1)*dim_loc;
    if(res >= lmin){ topo[0]=1; topo[1]=topo[2]; topo[2]=1; return test=1;}

    //on test si on peut intervertir la direction K vers I
    dim_loc = cells_tg/(nijk[1]*nijk[2]);
    res     = nijk[0]-(topo[2]-1)*dim_loc;
    if(res >= lmin){ topo[0]=topo[2]; topo[1]=1; topo[2]=1; return test=1;}
  }
  //
  // decoupe J
  //
  else if(topo[0]*topo[2] == 1)
  { //on test si on peut intervertir la direction J vers K
    dim_loc = cells_tg/(nijk[0]*nijk[1]);
    res     = nijk[2]-(topo[1]-1)*dim_loc;
    if(res >= lmin){ topo[2]=topo[1]; topo[0]=1; topo[1]=1; return test=1;}

    //on test si on peut intervertir la direction J vers I
    dim_loc = cells_tg/(nijk[1]*nijk[2]);
    res     = nijk[0]-(topo[1]-1)*dim_loc;
    if(res >= lmin){ topo[0]=topo[1]; topo[1]=1; topo[2]=1; return test=1;}
  }
  //
  // decoupe I
  //
  else if(topo[1]*topo[2] == 1)
  { //on test si on peut intervertir la direction I vers K
    dim_loc = cells_tg/(nijk[0]*nijk[1]);
    res     = nijk[2]-(topo[0]-1)*dim_loc;
    if(res >= lmin){ topo[2]=topo[0]; topo[0]=1; topo[1]=1; return test=1;}

    //on test si on peut intervertir la direction I vers J
    dim_loc = cells_tg/(nijk[0]*nijk[2]);
    res     = nijk[1]-(topo[0]-1)*dim_loc;
    if(res >= lmin){ topo[1]=topo[0]; topo[0]=1; topo[2]=1; return test=1;}
  }
  //
  // decoupe KJ
  //
  else if(topo[2] != 1)
  { //on test si on peut intervertir la direction J vers K et K vers J
    //J vers K
    dim_loc = cells_tg/(nijk[0]*nijk[1]/topo[2]);
    res     = nijk[2]-(topo[1]-1)*dim_loc;
    //printf("KKKKK %d %d %d %d %d \n", topo[0],topo[1],topo[2], res, dim_loc);
    if(res >= lmin && dim_loc >= lmin)
                   { //K vers J 
                    dim_loc = cells_tg/(nijk[0]*nijk[2]/topo[1]);
                    res     = nijk[1]-(topo[2]-1)*dim_loc;
                     //printf("JJJ %d %d %d %d %d \n", topo[0],topo[1],topo[2], res, dim_loc);
                    if(res >= lmin && dim_loc >= lmin ){ E_Int tmp=topo[2]; topo[2]=topo[1]; topo[1]=tmp; /*printf("cou1 \n");*/ return test=1; }
                   }
    
    //on test si on peut intervertir la direction J vers I et K vers J
    //J vers I
    dim_loc = cells_tg/(nijk[2]*nijk[1]/topo[0]);
    res     = nijk[0]-(topo[1]-1)*dim_loc;
    if(res >= lmin && dim_loc >= lmin ){ //K vers J 
                    dim_loc = cells_tg/(nijk[0]*nijk[2]/topo[1]);
                    res     = nijk[1]-(topo[2]-1)*dim_loc;
                    if(res >= lmin && dim_loc >= lmin ){ E_Int tmp=topo[0]; topo[0]=topo[1]; topo[1]= topo[2]; topo[2]=tmp;/*printf("cou2 \n");*/ return test=1;}
                   }

    //on test si on peut intervertir la direction J vers J et K vers I
    //J vers J
    dim_loc = cells_tg/(nijk[2]*nijk[0]);
    res     = nijk[1]-(topo[1]-1)*dim_loc;
    if(res >= lmin && dim_loc >= lmin){ //K vers I 
                    dim_loc = cells_tg/(nijk[1]*nijk[2]);
                    res     = nijk[0]-(topo[2]-1)*dim_loc;
                    if(res >= lmin && dim_loc >= lmin ){  topo[0]=topo[2]; topo[2]=1;  /*printf("cou3 \n");*/ return test=1;}
                   }
  }
  //
  // decoupe IJ
  //
  else if(topo[0] != 1)
  { //on test si on peut intervertir la direction J vers I et I vers J
    //J vers I
    dim_loc = cells_tg/(nijk[2]*nijk[1]);
    res     = nijk[0]-(topo[1]-1)*dim_loc;
    if(res >= lmin && dim_loc >= lmin){ //I vers J 
                    dim_loc = cells_tg/(nijk[0]*nijk[2]);
                    res     = nijk[1]-(topo[0]-1)*dim_loc;
                    if(res >= lmin && dim_loc >= lmin ){ E_Int tmp=topo[0]; topo[0]=topo[1]; topo[1]=tmp; topo[2]=1; return test=1;}
                   }
    
    //on test si on peut intervertir la direction I vers J et J vers K
    //I vers J
    dim_loc = cells_tg/(nijk[2]*nijk[0]);
    res     = nijk[1]-(topo[0]-1)*dim_loc;
    if(res >= lmin && dim_loc >= lmin ){ //J vers K 
                    dim_loc = cells_tg/(nijk[0]*nijk[1]);
                    res     = nijk[2]-(topo[1]-1)*dim_loc;
                    if(res >= lmin && dim_loc >= lmin ){ topo[2]=topo[1]; topo[1]=topo[0]; topo[0]=1; return test=1;}
                   }

    //on test si on peut intervertir la direction J vers J et I vers K
    //J vers J
    dim_loc = cells_tg/(nijk[2]*nijk[0]);
    res     = nijk[1]-(topo[1]-1)*dim_loc;
    if(res >= lmin && dim_loc >= lmin ){ //I vers K 
                    dim_loc = cells_tg/(nijk[1]*nijk[0]);
                    res     = nijk[2]-(topo[0]-1)*dim_loc;
                    if(res >= lmin && dim_loc >= lmin ){ topo[2]=topo[0]; topo[0]=1; return test=1;}
                   }
  }

 return test;

}

