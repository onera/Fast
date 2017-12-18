c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_ventijk( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro,
     &        param_int, param_real,
     &        ijkv_sdm,
     &        ind_sdm, ind_coe, ind_grad, 
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp,
     &        socket_topology, lok ,
     &        ti, tj, tk, vol, ti_df,tj_df,tk_df,
     &        venti, ventj, ventk, x, y ,z)

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
     & ithread, Nbre_socket, socket , compteur

      INTEGER_E  ijkv_sdm(3),ind_coe(6),ind_grad(6),ind_dm_zone(6),
     & ind_dm_omp(6), ind_dm_socket(6), ind_sdm(6), socket_topology(3),
     & param_int(0:*), lok(*)

      REAL_E ti(*),tj(*),tk(*),vol(*), venti(*),ventj(*),ventk(*),
     & x(*),y(*),z(*),ti_df(*),tj_df(*),tk_df(*)

      REAL_E param_real(0:*)
C Var loc 
      logical ksa, lerr
      INTEGER_E idir,nbdr_sdm,OMP_get_thread_num,icache,jcache,kcache,
     & Imax,Jmax,Kmax,itabi,ios,i,k, ind,lok_shap_loc,l,j,io,inc1,inc2,
     & size_max, size_loc,thread_parsock,inc11,inc22,l0,nitrun,nitcfg,
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

      nitrun   = -2
      nitcfg   = 1

#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"


            if(param_int(LALE).ge.1) then ! mise a jour Vent et tijk

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
     &                          venti, ventj, ventk, 
     &                          lok(1),lok(ipt_lok_sock),
     &                          lok(ipt_lok), omp_wait )



                call mjr_ale(ndo,nitcfg, ithread,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        x,y,z,ti,ti_df,tj,tj_df,tk,tk_df,vol,
     &                        venti, ventj, ventk)

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
     &                          venti, ventj, ventk, 
     &                          lok(1),lok(ipt_lok_sock),
     &                          lok(ipt_lok), omp_go )
              endif


#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      end
