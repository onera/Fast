c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine post( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, neq_grad,
     &        param_int, param_real, tke, enst, compteur,
     &        ijkv_sdm,
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
     & ithread, Nbre_socket, socket , neq_grad

      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6),ind_dm_omp(6),
     &  ind_dm_socket(6), socket_topology(3), param_int(0:*), lok(*)

      REAL_E rop(*), ti(*),tj(*),tk(*),vol(*), grad(*),tke,enst,compteur

      REAL_E param_real(0:*)

C Var loc 
      INTEGER_E nitrun
      INTEGER_E omp_mode, topo_omp(3), inddm_omp(6)
#include "FastS/HPC_LAYER/LOC_VAR_DECLARATION.for"
      REAL_E c1,c2,val

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !On force le mode 0 pour le Post. A supprimer quand mode 1 !stabiliser
      omp_mode = 0

      !coefficient pour calcul gradient ordre4 - c1=7./12 c2=1./12 ordre 4 - c1=0.5 c2=0 ordre 2
      c1 = 7./12
      c2 = 1./12

      tke      = 0.
      enst     = 0.
      compteur = 0 

      nitrun   = -2

#include "FastS/HPC_LAYER/SIZE_MIN.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"

           !Initilalisation systematique de grad
           val =0.
           call init_tab(ndo, val, param_int, param_int(NDIMDX),
     &                   neq_grad*3, ind_grad, grad )

#include "FastS/HPC_LAYER/SYNCHRO_WAIT.for"

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
     &                    tke  , enst , compteur,vol )

#include "FastS/HPC_LAYER/SYNCHRO_GO.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      end

