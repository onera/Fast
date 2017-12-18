c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine post( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, neq_grad,
     &        param_int, param_real, tke, enst, compteur,
     &        ijkv_sdm,
     &        ind_sdm, ind_coe, ind_grad, 
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp,
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

      INTEGER_E ndo, nidom, Nbre_thread_actif , mx_synchro, 
     & ithread, Nbre_socket, socket , neq_grad, compteur

      INTEGER_E  ijkv_sdm(3),ind_coe(6),ind_grad(6),ind_dm_zone(6),
     & ind_dm_omp(6), ind_dm_socket(6), ind_sdm(6), socket_topology(3),
     & param_int(0:*), lok(*)

      REAL_E rop(*), ti(*),tj(*),tk(*),vol(*), grad(*), tke,enst

      REAL_E param_real(0:*)
C Var loc 
      logical ksa, lerr
      INTEGER_E idir,nbdr_sdm,OMP_get_thread_num,icache,jcache,kcache,
     & Imax,Jmax,Kmax,itabi,ios,i,k, ind,lok_shap_loc,l,j,io,inc1,inc2,
     & size_max, size_loc,thread_parsock,inc11,inc22,l0,nitrun,
     & thread_parsock_actif,extended_range,
     & lok_shap_sock(4),thread_topology(3),
     & socket_pos(3), synchro_receive_sock(3),synchro_send_sock(3),
     & lok_shap(4), size_cache(3),
     & synchro_receive_th(3),synchro_send_th(3),ipt_lok_sock,
     & ithread_sock,ithread_io,ipt_lok,size_max_sock,neq_lok,taille,
     & ijkv_thread(3),kGbloc,jGbloc,iGbloc,ip,jp,kp,lth,
     & ibloc,jbloc,kbloc,ijkvloc(3),skip(3),shift(3),test(3),lwait,lgo,
     & size_thread(3),topo_s(3),thread_pos(3),thread_pos_tmp(3),sens(3),
     & ind_rhs(6),ind_mjr(6), size_target(3), cache(3)

      character*7 omp_init,omp_wait,omp_go,omp_wait_lu

      REAL_E c1,c2

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !coefficient pour calcul gradient ordre4 - c1=7./12 c2=1./12 ordre 4 - c1=0.5 c2=0 ordre 2
      c1 = 7./12
      c2 = 1./12

      tke      = 0.
      enst     = 0.
      compteur = 0.

      nitrun   = -2

#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"

           !Initilalisation systematique de grad
           call init_rhs(ndo, 1, param_int, param_int(NDIMDX),
     &                   neq_grad*3, ind_grad, grad )

            call synchro_omp_scater(param_int, ithread_io,
     &                          lth, sens,lgo,lwait,Nbre_socket,
     &                          Nbre_thread_actif,thread_parsock,
     &                          lok_shap_sock, lok_shap,neq_lok,
     &                          socket , socket_topology, socket_pos,
     &                          ithread, thread_topology,thread_pos_tmp,
     &                          synchro_receive_sock, synchro_send_sock,
     &                          synchro_receive_th  , synchro_send_th,
     &                          ibloc , jbloc , kbloc , ijkv_thread,
     &                          icache, jcache, kcache, ijkv_sdm,
     &                          ind_dm_omp,
     &                          grad, grad, grad, 
     &                          lok(1),lok(ipt_lok_sock),
     &                          lok(ipt_lok), omp_wait )


            call cp_gradu(ndo, ithread, neq_grad,
     &                    param_int, c1,c2,
     &                    ind_sdm,ind_mjr,
     &                    ind_dm_zone, ijkv_thread, ijkv_sdm,
     &                    synchro_send_sock, synchro_send_th,
     &                    synchro_receive_sock, synchro_receive_th,
     &                    ibloc , jbloc , kbloc ,
     &                    icache, jcache, kcache,
     &                    rop, grad,ti,tj,tk,vol)

            call cptaylor(ndo, neq_grad,
     &                    param_int, ind_mjr,
     &                    rop  , grad, 
     &                    tke  , enst , compteur )



            call synchro_omp_scater(param_int, ithread_io,
     &                          lth, sens,lgo,lwait,Nbre_socket,
     &                          Nbre_thread_actif,thread_parsock,
     &                          lok_shap_sock, lok_shap,neq_lok,
     &                          socket , socket_topology, socket_pos,
     &                          ithread, thread_topology,thread_pos_tmp,
     &                          synchro_receive_sock, synchro_send_sock,
     &                          synchro_receive_th  , synchro_send_th,
     &                          ibloc , jbloc , kbloc , ijkv_thread,
     &                          icache, jcache, kcache, ijkv_sdm,
     &                          ind_dm_omp,
     &                          grad, grad, grad, 
     &                          lok(1),lok(ipt_lok_sock),
     &                          lok(ipt_lok), omp_go )


#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      end
