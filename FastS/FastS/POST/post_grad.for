c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine post_grad( ndo, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, neq_grad, order,
     &        param_int, param_real, 
     &        ijkv_sdm,
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp, topo_s,
     &        socket_topology, lok ,
     &        rop , ti, tj, tk, vol, grad)

c***********************************************************************
c_U   USER : TERRACOL
c
c     ACT
c_A    Appel du calcul des flux explicites
c
c     VAL
c_V    gaz parfait monoespece
c_V    processeur domaine
c_V    steady/unsteady
c
c     INP
c_I    tijk     : vecteur param_int( IO_THREAD)rmale aux facettes des mailles
c_I    vent     : vitesses d'entrainement aux facettes preced.
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    grad    : increment de la solution
c
c***********************************************************************
      implicit none

      INTEGER_E ndo, Nbre_thread_actif , mx_synchro, 
     & ithread, Nbre_socket, socket , neq_grad, order

      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6), topo_s(3),
     & ind_dm_omp(6), ind_dm_socket(6), socket_topology(3),
     & param_int(0:*), lok(*)

      REAL_E rop(*), ti(*),tj(*),tk(*),vol(*), grad(*)

      REAL_E param_real(0:*)
C Var loc 
      INTEGER_E nitrun,V1,V2,V3, l,lij,j,k,ind_loop(6),ltij,lt,lvo
#include "../FastC/FastC/HPC_LAYER/LOC_VAR_DECLARATION.for"

      REAL_E c1,c2,volinv, val

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !coeficient pour calcul gradient ordre4 !c1=0.5 c2 =0 ordre 2
      if(order.eq.4) then
        c1 = 7./6.
        c2 = 1./6.
      ! ordre 2 obligatoire si calcul sur une rangee fictive
      else
        c1 = 1.
        c2 = 0.
      endif

      nitrun =-2

      V1 = 0
      V2 =   param_int(NDIMDX)
      V3 = 2*param_int(NDIMDX)

#include "../FastC/FastC/HPC_LAYER/SIZE_MIN.for"
#include "../FastC/FastC/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
      if(c1.eq.1.) extended_range = 1
#include "../FastC/FastC/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "../FastC/FastC/HPC_LAYER/INDICE_RANGE.for"

           !Initilalisation systematique de grad
           val =0.
           call init_tab(ndo, val, param_int, param_int(NDIMDX),
     &                   neq_grad, ind_rhs, grad )

#include "../FastC/FastC/HPC_LAYER/SYNCHRO_WAIT.for"

            if(param_int(ITYPZONE).eq.0) then

               call cp_gradvar_3dfull(ndo, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        rop, grad,ti,tj,tk)
            
                ind_loop(:) = ind_mjr(:)
#include        "FastS/Compute/loop_begin.for"
                   volinv=0.5/vol(lvo)
                   grad(l+ V1) = grad(l+ V1)*volinv
                   grad(l+ V2) = grad(l+ V2)*volinv
                   grad(l+ V3) = grad(l+ V3)*volinv
#include         "FastS/Compute/loop_end.for"

            elseif(param_int(ITYPZONE).eq.1) then

               call cp_gradvar_3dhomo(ndo, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        rop, grad,ti,tj,tk)
              
                ind_loop(:) = ind_mjr(:)
#include        "FastS/Compute/loop_begin.for"
                   volinv=0.5/vol(lvo)
                   grad(l+ V1) = grad(l+ V1)*volinv
                   grad(l+ V2) = grad(l+ V2)*volinv
                   grad(l+ V3) = grad(l+ V3)*volinv
#include         "FastS/Compute/loop_end.for"

            elseif(param_int(ITYPZONE).eq.2) then

              
               call cp_gradvar_3dcart(ndo, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        rop, grad,ti,tj,tk)

                 volinv=0.5/vol(1)
                 ind_loop(:) = ind_mjr(:)
#include        "FastS/Compute/loop_begin.for"
                   grad(l+ V1) = grad(l+ V1)*volinv
                   grad(l+ V2) = grad(l+ V2)*volinv
                   grad(l+ V3) = grad(l+ V3)*volinv
#include         "FastS/Compute/loop_end.for"


            else
               call cp_gradvar_2d(ndo, ithread, neq_grad,
     &                        param_int, c1,c2,
     &                        ind_sdm,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        rop, grad,ti,tj,tk)

                 ind_loop(:) = ind_mjr(:)
#include        "FastS/Compute/loop_begin.for"
                   volinv=0.5/vol(lvo)
                   grad(l+ V1) = grad(l+ V1)*volinv
                   grad(l+ V2) = grad(l+ V2)*volinv
#include         "FastS/Compute/loop_end.for"

            endif

            call extrap(ndo,param_int,c1,ind_sdm,ind_dm_zone,grad(1+V1))
            call extrap(ndo,param_int,c1,ind_sdm,ind_dm_zone,grad(1+V2))
            if(param_int(ITYPZONE).ne.3) 
     &      call extrap(ndo,param_int,c1,ind_sdm,ind_dm_zone,grad(1+V3))
               
#include "../FastC/FastC/HPC_LAYER/SYNCHRO_GO.for"
#include "../FastC/FastC/HPC_LAYER/LOOP_CACHE_END.for"
#include "../FastC/FastC/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      end
