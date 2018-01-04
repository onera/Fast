c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine post_qprime_1rot( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, 
     &        order, dim_grad, vr1, vr2,
     &        param_int, param_real,
     &        ijkv_sdm,
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp,
     &        socket_topology, lok ,
     &        rop , rop_m1, ti, tj, tk, vol, Q, rot1)

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

      INTEGER_E ndo, nidom, Nbre_thread_actif , mx_synchro, order,
     & ithread, Nbre_socket, socket , dim_grad, vr1, vr2

      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6),ind_dm_omp(6),
     &  ind_dm_socket(6), socket_topology(3), param_int(0:*), lok(*)

      REAL_E rop(*), rop_m1(*), ti(*),tj(*),tk(*),vol(*), Q(*), rot1(*)

      REAL_E param_real(0:*)

C Var loc 
      INTEGER_E nitrun
#include "FastS/HPC_LAYER/LOC_VAR_DECLARATION.for"
      REAL_E c1,c2

      REAL_E dvardc(dim_grad*3*3)

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"

      !coeficient pour calcul gradient ordre2 !
      !coeficient pour calcul gradient ordre4 !c1=0.5 c2 =0 ordre 2
      if(order.eq.4) then
        c1 = 7./6.
        c2 = 1./6.
      ! ordre 2 obligatoire si calcul sur une rangee fictive
      else
        c1 = 1.0
        c2 = 0.
      endif

      nitrun = -2

#include "FastS/HPC_LAYER/SIZE_MIN.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
      if(c1.eq.1.) extended_range = 1
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"

            if(param_int(ITYPZONE).eq.0) then

               call cp_qprime_1rot_3dfull(ndo,ithread, dim_grad,vr1,vr2,
     &                        param_int, param_real, c1,c2,
     &                        ind_sdm ,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        dvardc,
     &                        rop, rop_m1, Q,rot1,ti,tj,tk, vol)
            
            elseif(param_int(ITYPZONE).eq.1) then

               call cp_qprime_1rot_3dhomo(ndo,ithread, dim_grad,vr1,vr2,
     &                        param_int, param_real, c1,c2, 
     &                        ind_sdm ,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        dvardc,
     &                        rop, rop_m1, Q, rot1, ti,tj,tk, vol)

            elseif(param_int(ITYPZONE).eq.2) then

              
               call cp_qprime_1rot_3dcart(ndo,ithread, dim_grad,vr1,vr2,
     &                        param_int, param_real, c1,c2,
     &                        ind_sdm ,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        dvardc,
     &                        rop, rop_m1, Q, rot1, ti,tj,tk, vol)

            else
               call cp_qprime_2d(ndo, ithread,  dim_grad,vr1,vr2,
     &                        param_int, param_real, c1,c2,
     &                        ind_sdm ,
     &                        ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        dvardc,
     &                        rop, rop_m1, Q, rot1, ti,tj,tk, vol)
            endif


            call extrap(ndo, param_int, c1, ind_sdm, ind_dm_zone, Q   )
            call extrap(ndo, param_int, c1, ind_sdm, ind_dm_zone, rot1)

#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      end
