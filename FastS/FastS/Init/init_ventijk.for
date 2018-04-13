c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_ventijk( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, omp_mode,
     &        param_int, param_real,
     &        ijkv_sdm,
     &        ind_dm_zone, ind_dm_socket,
     &        socket_topology, lok , topo_omp, inddm_omp,
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
     & ithread, Nbre_socket, socket , compteur,omp_mode

      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6),
     & ind_dm_omp(6), ind_dm_socket(6), socket_topology(3),
     & param_int(0:*), lok(*), topo_omp(3), inddm_omp(6)

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

      if(param_int(LALE).ge.1) then ! mise a jour Vent et tijk


        !coefficient pour calcul gradient ordre4 - c1=7./12 c2=1./12 ordre 4 - c1=0.5 c2=0 ordre 2
        c1 = 7./12
        c2 = 1./12

        nitrun   = -2
        nitcfg   = 1
#include "FastS/HPC_LAYER/SIZE_MIN.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"
#include "FastS/HPC_LAYER/SYNCHRO_WAIT.for"

                call mjr_ale(ndo,nitcfg, ithread,
     &                        param_int, param_real,
     &                        ind_dm_zone, ind_sdm,ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th,
     &                        synchro_receive_sock, synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        x,y,z,ti,ti_df,tj,tj_df,tk,tk_df,vol,
     &                        venti, ventj, ventk)

#include "FastS/HPC_LAYER/SYNCHRO_GO.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      endif

      end
