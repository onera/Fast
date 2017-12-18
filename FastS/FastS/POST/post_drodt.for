c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine post_drodt( ndo, nidom, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, flag, dim_grad,
     &        param_int, param_real,
     &        ijkv_sdm,
     &        ind_sdm, ind_coe, ind_grad, 
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp,
     &        socket_topology, lok ,
     &        rop , rop_n, rop_n1, drodt)
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

      INTEGER_E ndo, nidom, Nbre_thread_actif , mx_synchro, flag,
     & ithread, Nbre_socket, socket , dim_grad

      INTEGER_E  ijkv_sdm(3),ind_coe(6),ind_grad(6),ind_dm_zone(6),
     & ind_dm_omp(6), ind_dm_socket(6), ind_sdm(6), socket_topology(3),
     & param_int(0:*), lok(*)

      REAL_E rop(*), rop_n(*), rop_n1(*), drodt(*)

      REAL_E param_real(0:*)
C Var loc 
      logical ksa, lerr
      INTEGER_E idir,nbdr_sdm,OMP_get_thread_num,icache,jcache,kcache,
     & Imax,Jmax,Kmax,itabi,ios,i,k, ind,lok_shap_loc,l,j,io,inc1,inc2,
     & size_max, size_loc,thread_parsock,inc11,inc22,l0,nitrun,
     & thread_parsock_actif,lij,extended_range,
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

      REAL_E c1,c2,c3, dtinv

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"

      !coeficient pour calcul gradient ordre2 !
      dtinv = 1./param_real(DTC)

      c1 =  1.5*dtinv
      c2 = -2.*dtinv
      c3 =  0.5*dtinv

      nitrun = -2

#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for"
#include "FastS/HPC_LAYER/LOOP_CACHE_BEGIN.for"
#include "FastS/HPC_LAYER/INDICE_RANGE.for"

#ifndef E_SCALAR_COMPUTER
          incmax=param_int(NIJK)*param_int(NIJK+1)*param_int(NIJK+4)
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
          do 100 l=1,param_int(NDIMDX)
#else
        do 100 k = ind_coe(5), ind_coe(6)
        do 100 j = ind_coe(3), ind_coe(4)
           lij  = inddm(ind_coe(1) , j, k)
CDIR$ IVDEP
        do 100 l = lij, lij +  ind_coe(2) - ind_coe(1)
#endif
            drodt(l)= c1*rop(l) + c2*rop_n(l)+ c3*rop_n1(l)

  100 continue

#include "FastS/HPC_LAYER/LOOP_CACHE_END.for"
#include "FastS/HPC_LAYER/WORK_DISTRIBUTION_END.for"

      end
