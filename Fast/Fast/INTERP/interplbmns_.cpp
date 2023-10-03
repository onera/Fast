/*
    Copyright 2013-2017 Onera.

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
# include <iostream>
# include <string.h>
# include "../FastS/FastS/fastS.h"
# include "fast.h"
# include "../FastASLBM/FastASLBM/fastLBM.h"
# include "param_solver.h"
# include <string.h>
# include <omp.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Idem: in place + from zone + tc compact au niveau base
//=============================================================================
PyObject* K_FAST::interplbmns_(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *work;
  PyObject *pyParam_int, *pyParam_real;
  E_Int loc, nstep, nitrun, nitmax;
  E_Int NoTransfert, process;

  if (!PYPARSETUPLE_(args, OOOO_ O_ IIII_ I_, 
                    &zonesR, &zonesD, &pyParam_int, &pyParam_real,&work, &nitmax, &nitrun, &nstep,  &NoTransfert, &process))
  {
      return NULL;
  }

  PyObject* dtlocArray  = PyDict_GetItemString(work,"dtloc"); FldArrayI* dtloc;
  K_NUMPY::getFromNumpyArray(dtlocArray, dtloc, true); E_Int* iptdtloc  = dtloc->begin();
  E_Int nssiter = iptdtloc[0];
  E_Int omp_mode = iptdtloc[8];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;


  //Nombre de zones doneuses et receveuses
  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //// Recuperation du tableau param_int de l'arbre t et des veceurs rop et roptmp
  // Transferts se font sur Density et Density_P1
  E_Int** param_int    = new E_Int*[nidomR];
  E_Float** param_real = new E_Float*[nidomR*3];
  E_Float** iptrom_CL  = param_real + nidomR; 
  E_Float** iptrom     = iptrom_CL  + nidomR;

  vector<PyArrayObject*> hook;
  // On parcourt les zones R et on recupere les tableaux/vecteurs
  for (E_Int nd = 0; nd < nidomR; nd++)
     {
       PyObject* zone = PyList_GetItem(zonesR, nd);

       PyObject* numerics= K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
       PyObject*      o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"  );
       param_int[nd]     = K_PYTREE::getValueAI(o, hook);

	              o = K_PYTREE::getNodeFromName1(numerics, "Parameter_real" );
       param_real[nd]   = K_PYTREE::getValueAF(o, hook);

                      o = K_PYTREE::getNodeFromName1(zone     , "FlowSolution#Centers");

       PyObject*      t = K_PYTREE::getNodeFromName1( o       , "Density"     );
       iptrom[nd]       = K_PYTREE::getValueAF(t, hook);

       if( (param_int[nd][ITYPCP]<=1) || (param_int[nd][ITYPCP]==2 && nstep!=2) )
         {  PyObject*      t = K_PYTREE::getNodeFromName1( o       , "Density_P1"     );
            iptrom_CL[nd]    = K_PYTREE::getValueAF(t, hook);
         }
       else {  iptrom_CL[nd]    = iptrom[nd]; }
     }

//#pragma omp parallel default(shared) num_threads(1)//private(cycle)
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num() +1;
    E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
    E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
    if( Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

    E_Int Nbre_thread_actif_loc, ithread_loc;
    if( omp_mode == 1) { Nbre_thread_actif_loc = 1;                 ithread_loc = 1;}
    else               { Nbre_thread_actif_loc = Nbre_thread_actif; ithread_loc = ithread;}

  E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
  E_Int socket          = (ithread-1)/thread_parsock +1;
  E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

  E_Float rhs_end=0;

  E_Int nitcfg_loc= 2;

  E_Int nbtask = ipt_omp[nitcfg_loc-1]; 
  E_Int ptiter = ipt_omp[nssiter+ nitcfg_loc-1];

  //calcul du sous domaine a traiter par le thread
  for (E_Int ntask = 0; ntask < nbtask; ntask++)
     {
       E_Int pttask = ptiter + ntask*(6+Nbre_thread_actif*7);
       E_Int nd = ipt_omp[ pttask ];

       E_Int ndo   = nd;

       E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp;

       ithread_loc           = ipt_omp[ pttask + 2 + ithread -1 ] +1 ;
       E_Int nd_subzone      = ipt_omp[ pttask + 1 ];
       Nbre_thread_actif_loc = ipt_omp[ pttask + 2 + Nbre_thread_actif ];
       ipt_topo_omp          = ipt_omp + pttask + 3 + Nbre_thread_actif ;
       ipt_inddm_omp         = ipt_omp + pttask + 2 + Nbre_thread_actif +4 + (ithread_loc-1)*6;

       if (ithread_loc == -1) {continue;}

       if (param_int[nd][ ITYPCP ] <= 1) {continue;} //pas d'interp/stck a gerer si NS implicit

       E_Int pt_interp = param_int[nd][PT_INTERP];
       E_Int nrac = param_int[nd][pt_interp];
       for (E_Int rac=0; rac < nrac; rac++) // Boucle sur les differents raccords
         {
            E_Int pt_racInt = param_int[nd][ pt_interp + rac +1        ]; //pointer vers info Int du raccord (type, range)
            E_Int pt_racReal= param_int[nd][ pt_interp + rac +1 +nrac  ]; //pointer vers stockage float
             
            E_Int typ_rac =  param_int[nd][pt_racInt];
            if(typ_rac == 1)
             {
               E_Int* range  =  param_int[nd] + pt_racInt + 3;
               E_Int  size   =  (range[1]-range[0]+1)* (range[3]-range[2]+1)* (range[5]-range[4]+1);
             
               E_Int range_th[6];

               range_th[0]= 1; range_th[1]= 0; range_th[2]= 1;range_th[3]= 0;range_th[4]= 1;range_th[5]= 0;

               //intersection entre plage thread et fenetre stockage
               if(   ipt_inddm_omp[0]<= range[1] &&  ipt_inddm_omp[1] >= range[0] 
                  && ipt_inddm_omp[2]<= range[3] &&  ipt_inddm_omp[3] >= range[2] 
                  && ipt_inddm_omp[4]<= range[5] &&  ipt_inddm_omp[5] >= range[4] )
                 {
                   range_th[0]=max( range[0] , ipt_inddm_omp[0]);
                   range_th[1]=min( range[1] , ipt_inddm_omp[1]);
                   range_th[2]=max( range[2] , ipt_inddm_omp[2]);
                   range_th[3]=min( range[3] , ipt_inddm_omp[3]);
                   range_th[4]=max( range[4] , ipt_inddm_omp[4]);
                   range_th[5]=min( range[5] , ipt_inddm_omp[5]);
                 }

                 if(nd==1) for (E_Int i=0; i < 6; i++) { printf("rg %d %d %d  th= %d , idir= %d \n", range_th[i], range[i] , ipt_inddm_omp[i], ithread_loc, i+1);}

               E_Float* iptS = NULL; E_Float* iptPsiG = NULL; E_Float* iptro_CL = NULL;
               if(nstep <=2)
                 {
                  E_Int sens=0; //(ro --> stk)
                  E_Int pos_p1 = param_int[nd][pt_racInt +1];
                  E_Int   neq  = param_int[nd][pt_racInt +2];
                  pos_p1+=1; if(pos_p1==3){pos_p1=0;}

                  E_Float* stk  = param_real[nd]+  pt_racReal + pos_p1*neq*size;

                  printf("COPY: %d rac= %d , nd= %d %d , pos= %d , th= %d , nstep= %d \n ",  pt_racReal + pos_p1*neq*size, rac, nd, neq, pos_p1, ithread_loc, nstep);
                  //sauvegarde de macro (instant N) dans stk 
                  //                                                              macro,         distribution,   Sij,   , autre terme...
                  if(nstep==1){copy_values_( nd, neq, param_int[nd], range, range_th, size,  stk, iptrom[nd], iptro_CL,  iptS, iptPsiG, sens, typ_rac);}

                  stk  = param_real[nd]+  pt_racReal;
                  //                                                             macro,         distribution,   Sij,   , autre terme...
                  interp_lbm_dtloc_( neq, param_int[nd], range, range_th, size,  stk, iptrom_CL[nd], iptro_CL,  iptS, iptPsiG, pos_p1, typ_rac, nstep);
                 }
                else
                 {
                  E_Int sens=1; //(stk--> ro_p1)
                  E_Int pos_p1 = param_int[nd][pt_racInt +1] + 1;
                  if(pos_p1==3){pos_p1=0;}

                  E_Int   neq  = param_int[nd][pt_racInt +2];
           
                  printf("Recup: %d rac= %d , nd= %d %d , pos= %d , th= %d , nstep= %d \n ",  pt_racReal + pos_p1*neq*size, rac, nd, neq, pos_p1, ithread_loc, nstep);
                  E_Float* stk  = param_real[nd]+  pt_racReal + pos_p1*neq*size;
     
                  //                                                             macro,         distribution,   Sij,   , autre terme...
                  copy_values_( nd, neq, param_int[nd], range, range_th, size,  stk, iptrom[nd], iptro_CL,  iptS, iptPsiG, sens, typ_rac);
                 }

             }//typ_rac
         }//loop rac


     }//loop task
  }//omp fin

 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
