c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_ventijk( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro,
     &        param_int, param_real,
     &        ijkv_sdm,
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

      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6),
     & ind_dm_omp(6), ind_dm_socket(6), socket_topology(3),
     & param_int(0:*), lok(*)

      REAL_E ti(*),tj(*),tk(*),vol(*), venti(*),ventj(*),ventk(*),
     & x(*),y(*),z(*),ti_df(*),tj_df(*),tk_df(*)

      REAL_E param_real(0:*)
C Var loc 
      INTEGER_E nitrun,nitcfg
#include "FastS/HPC_LAYER/LOC_VAR_DECLARATION.for"

      REAL_E c1,c2

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !coefficient pour calcul gradient ordre4 - c1=7./12 c2=1./12 ordre 4 - c1=0.5 c2=0 ordre 2
      c1 = 7./12
      c2 = 1./12

      nitrun   = -2
      nitcfg   = 1
#include "FastS/HPC_LAYER/SIZE_MIN.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"


            if(param_int(LALE).ge.1) then ! mise a jour Vent et tijk

            call synchro_omp_scater(param_int, ithread,
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

            call synchro_omp_scater(param_int, ithread,
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
