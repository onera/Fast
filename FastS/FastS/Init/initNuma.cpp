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
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
//# include <numa.h>
//# include </usr/src/kernels/2.6.32-573.12.1.el6.x86_64/include/linux/numa.h>

# include "fastS.h"
//# include "converter.h"
# include "kcore.h"
# include "param_solver.h"

using namespace std;
using namespace K_FLD;

extern "C"
{
  E_Int numa_available();  E_Int numa_move_pages( int i, unsigned long j, void** pt, const int* node,  int* stat, int err);
}
//=============================================================================
/* Init data in parallel openmp to improve data placement on numa machine */
//=============================================================================
PyObject* K_FASTS::initNuma(PyObject* self, PyObject* args)
{
  PyObject* sourceArray;  PyObject* targetArray;  PyObject* param_int; E_Int ivar;  E_Int thread_numa;

#ifdef E_DOUBLEINT 
  if (!PyArg_ParseTuple(args, "OOOll", &sourceArray, &targetArray, &param_int, &ivar, &thread_numa ))
#else
  if (!PyArg_ParseTuple(args, "OOOii", &sourceArray, &targetArray, &param_int, &ivar, &thread_numa ))
#endif
  {
    return NULL;
  }



  FldArrayF* source; FldArrayF* target;
  K_NUMPY::getFromNumpyArray(sourceArray, source, true); E_Float* iptsource = source->begin();
  K_NUMPY::getFromNumpyArray(targetArray, target, true); E_Float* ipttarget = target->begin();


  vector<PyArrayObject*> hook;
  E_Int* ipt_param_int  = K_PYTREE::getValueAI(param_int, hook);

  E_Int ndo    =1;
  E_Int max_thread = 1;
#ifdef _OPENMP
  max_thread = omp_get_max_threads(); // !nombre de thread maximal dans le calcul
#endif

  FldArrayI ind_dm(6); E_Int* ipt_ind_dm  = ind_dm.begin();

  
  ipt_ind_dm[0] = 1;
  ipt_ind_dm[2] = 1;
  ipt_ind_dm[4] = 1;
  ipt_ind_dm[1] = ipt_param_int[ IJKV   ]; 
  ipt_ind_dm[3] = ipt_param_int[ IJKV +1];
  ipt_ind_dm[5] = ipt_param_int[ IJKV +2];

  FldArrayI thread_topology(3*max_thread);
  FldArrayI ind_dm_socket(  6*max_thread);  
  FldArrayI ind_dm_thread(  6*max_thread);  
  FldArrayI      ind_loop(  6*max_thread);  

//     E_Int shiftvar = ivar* ipt_param_int[NIJK]*ipt_param_int[NIJK+1]*ipt_param_int[NIJK+2];
//  printf("ni =%d %d %d %d \n",ipt_param_int[0],ipt_param_int[1],ipt_param_int[2], shiftvar);
//  printf("ni =%d %d %d %d  \n", ipt_ind_dm[1],ipt_ind_dm[3],ipt_ind_dm[5],thread_numa);

#pragma omp parallel default(shared)
  {


#ifdef _OPENMP 
       E_Int  ithread           = omp_get_thread_num() +1;
       E_Int  Nbre_thread_actif = omp_get_num_threads(); // !nombre de thread actif dans cette zone
#else
       E_Int  ithread           = 1;
       E_Int  Nbre_thread_actif = 1; // !nombre de thread actif dans cette zone
#endif
     E_Int Nbre_socket   = NBR_SOCKET;     // nombre socket sur le noeud a memoire partagee
     if( Nbre_thread_actif < Nbre_socket ) Nbre_socket = 1;

     E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
     E_Int socket          = (ithread-1)/thread_parsock +1;
     E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

 
     E_Int shiftvar = ivar* (ipt_param_int[NIJK]*ipt_param_int[NIJK+1]*ipt_param_int[NIJK+2] + ipt_param_int[SHIFTVAR]) ;

  if(thread_numa == 0)
  { 
     E_Int* ipt_thread_topology   = thread_topology.begin()   + 3*(ithread-1);
     E_Int* ipt_ind_dm_socket     = ind_dm_socket.begin()     + 6*(ithread-1);
     E_Int* ipt_ind_dm_thread     = ind_dm_thread.begin()     + 6*(ithread-1);
     E_Int* ipt_ind_loop          = ind_loop.begin()          + 6*(ithread-1);

     //decoupage socket
     E_Int lmin = 10;
     if (ipt_param_int[ ITYPCP ] == 2) lmin = 4;

     indice_boucle_lu_(ndo, socket , Nbre_socket, lmin, ipt_ind_dm, 
                       ipt_thread_topology, ipt_ind_dm_socket );

     //decoupage sur les threads du socket
     indice_boucle_lu_(ndo, ithread_sock, thread_parsock, lmin, ipt_ind_dm_socket, 
                       ipt_thread_topology, ipt_ind_dm_thread );


     E_Int ific   = ipt_param_int[ NIJK +3];
     if(ipt_ind_dm[1]==1) ific = 0;   // prise en compte moyenne tmy dir I homogene
     E_Int jfic   = ipt_param_int[ NIJK +3];
     if(ipt_ind_dm[3]==1) jfic = 0;   // prise en compte moyenne tmy dir J homogene
     E_Int kfic   = ipt_param_int[ NIJK +4];
     if(ipt_ind_dm[5]==1) kfic = 0;   // prise en compte moyenne tmy dir K homogene

     E_Int iloop1 = ipt_ind_dm_thread[0];
     E_Int jloop1 = ipt_ind_dm_thread[2];
     E_Int kloop1 = ipt_ind_dm_thread[4];
     if( iloop1 == 1) iloop1 = iloop1 - ific;
     if( jloop1 == 1) jloop1 = jloop1 - jfic;
     if( kloop1 == 1) kloop1 = kloop1 - kfic;

     E_Int iloop2 = ipt_ind_dm_thread[1];
     E_Int jloop2 = ipt_ind_dm_thread[3];
     E_Int kloop2 = ipt_ind_dm_thread[5];
     if( iloop2 == ipt_ind_dm[1]) iloop2 = iloop2 +  ific;
     if( jloop2 == ipt_ind_dm[3]) jloop2 = jloop2 +  jfic ;
     if( kloop2 == ipt_ind_dm[5]) kloop2 = kloop2 +  kfic;

     E_Int ni_loc = ipt_param_int[ NIJK   ];
     E_Int nj_loc = ipt_param_int[ NIJK +1];

//  printf("ific =%d %d  %d %d %d  % d %d %d \n", ipt_thread_topology[0], ipt_thread_topology[1], ipt_thread_topology[2],jloop2, ific,kfic, ni_loc,nj_loc );
  //if(ithread==1) printf("ptr = %p \n", ipttarget);
    
    // Distribution de la sous-zone sur les threads
    ipt_ind_loop[0]= iloop1; ipt_ind_loop[1]= iloop2;
    ipt_ind_loop[2]= jloop1; ipt_ind_loop[3]= jloop2;
    ipt_ind_loop[4]= kloop1; ipt_ind_loop[5]= kloop2;

    //copynuma_(ipt_ind_loop, ni_loc, nj_loc, shiftvar, ific,jfic,kfic, ipttarget, iptsource);

      for ( E_Int k = kloop1; k <= kloop2; k++)
        {
         for ( E_Int j = jloop1; j <= jloop2; j++)
           {
            for ( E_Int i = iloop1; i <= iloop2; i++)
             {

               E_Int lsrc = (i+ific-1) + (j+jfic-1)*ni_loc +(k+kfic-1)*ni_loc*nj_loc;
               E_Int l    = lsrc +shiftvar;

               ipttarget[l ] = iptsource[lsrc];
             }
           }
        }
      #pragma omp barrier
/*   if(ithread >=1){ E_Int ibloc = (kloop2-kloop1+1)*(jloop2-jloop1+1)*(iloop2-iloop1+1)/512;

                    int status[1];
                    int node[1]; node[0]=1;
                    void *pa;
                    unsigned long a;
                    // round p down to the nearest page boundary

                   a  = (unsigned long) ipttarget;
                   a   = a - (a % ((unsigned long) getpagesize()));
                   pa  = (void *) a;
                   //printf("ptro   = %p %d %d %d %d %d\n", pa, getpagesize(), ibloc, numa_available(), ithread, MPOL_MF_MOVE );
                   //printf("ptro   = %p %d %d %d %d \n", pa, getpagesize(), ibloc, numa_available(), ithread );
                   for ( E_Int k = 0; k < ibloc; k++)
                            {
                             a  = (unsigned long) ipttarget+ k*512*8;
                             unsigned long b   = a - (a % ((unsigned long) getpagesize()));
                             pa  = (void *) b;
                             //if(ithread==17) printf("a,b=   = %d %d %d  \n", a, b, k);
                             //numa_move_pages( getpid(),1,&pa,node,status,0) ;
                             if (numa_move_pages(getpid(),1,&pa,node,status,0) != 0) { printf("Problem in calling move_pages()\n");}
                             //if(ithread==17) printf("ptro   = %p %d %d %d \n", pa, status[0], k, getpid());
                            }
                   } 
 */
  }
  else //omp_mode=1           
  {

    FldArrayI ind_dm(6); E_Int* inddm  = ind_dm.begin();

    E_Int shift_omp = ipt_param_int[ PT_OMP ];

    E_Int Nbre_thread_actif_loc = ipt_param_int[ shift_omp  + Nbre_thread_actif ];
    E_Int ithread_loc           = ipt_param_int[ shift_omp  +  ithread -1       ] +1 ;


    E_Int* inddm_omp = ipt_param_int + shift_omp +  Nbre_thread_actif + 4 + (ithread_loc-1)*6;

    inddm[0] = inddm_omp[0];
    inddm[1] = inddm_omp[1];
    inddm[2] = inddm_omp[2];
    inddm[3] = inddm_omp[3];
    inddm[4] = inddm_omp[4];
    inddm[5] = inddm_omp[5];

    E_Int iloop1 = K_FUNC::E_max(1, inddm[0] );              //prise en compte moyenne tmy dir I homogene
    E_Int jloop1 = K_FUNC::E_max(1, inddm[2] );              //prise en compte moyenne tmy dir I homogene
    E_Int kloop1 = K_FUNC::E_max(1, inddm[4] );              //prise en compte moyenne tmy dir I homogene
    E_Int iloop2 = K_FUNC::E_min(ipt_ind_dm[1], inddm[1] );  //prise en compte moyenne tmy dir I homogene
    E_Int jloop2 = K_FUNC::E_min(ipt_ind_dm[3], inddm[3] );  //prise en compte moyenne tmy dir I homogene
    E_Int kloop2 = K_FUNC::E_min(ipt_ind_dm[5], inddm[5] );  //prise en compte moyenne tmy dir I homogene

     E_Int ific   = ipt_param_int[ NIJK +3];
     if(ipt_ind_dm[1]==1) ific = 0;   // prise en compte moyenne tmy dir I homogene
     E_Int jfic   = ipt_param_int[ NIJK +3];
     if(ipt_ind_dm[3]==1) jfic = 0;   // prise en compte moyenne tmy dir J homogene
     E_Int kfic   = ipt_param_int[ NIJK +4];
     if(ipt_ind_dm[5]==1) kfic = 0;   // prise en compte moyenne tmy dir K homogene

     if(iloop1 ==1) iloop1 = iloop1 -ific;
     if(jloop1 ==1) jloop1 = jloop1 -ific;
     if(kloop1 ==1) kloop1 = kloop1 -kfic;
     if(iloop2 ==ipt_ind_dm[1]) iloop2 = iloop2 +ific;
     if(jloop2 ==ipt_ind_dm[3]) jloop2 = jloop2 +ific;
     if(kloop2 ==ipt_ind_dm[5]) kloop2 = kloop2 +kfic;

     E_Int ni_loc = ipt_param_int[ NIJK   ];
     E_Int nj_loc = ipt_param_int[ NIJK +1];
  
     //printf("ific =%d %d  %d %d %d  % d %d %d \n", iloop1, iloop2, jloop1, jloop2,  kloop1, kloop2, ific,kfic  );
     if(ithread_loc != -1)
     {
      for ( E_Int k = kloop1; k <= kloop2; k++)
        {
         for ( E_Int j = jloop1; j <= jloop2; j++)
           {
            for ( E_Int i = iloop1; i <= iloop2; i++)
             {

               E_Int lsrc = (i+ific-1) + (j+jfic-1)*ni_loc +(k+kfic-1)*ni_loc*nj_loc;
               E_Int l    = lsrc +shiftvar;

               ipttarget[l ] = iptsource[lsrc];
             }
           }
        }
     }
  }// mode_omp

  }// omp


  RELEASESHAREDN( sourceArray  , source  );
  RELEASESHAREDN( targetArray  , target  );
  RELEASEHOOK(hook)

  Py_INCREF(Py_None);
  return Py_None;

}
