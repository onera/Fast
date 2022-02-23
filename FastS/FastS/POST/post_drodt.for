c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine post_drodt( ndo, Nbre_thread_actif, 
     &        ithread, Nbre_socket, socket, mx_synchro, flag, dim_grad,
     &        param_int, param_real,
     &        ijkv_sdm,
     &        ind_dm_zone, ind_dm_socket, ind_dm_omp, topo_s,
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

      INTEGER_E ndo, Nbre_thread_actif , mx_synchro, flag,
     & ithread, Nbre_socket, socket , dim_grad

      INTEGER_E  ijkv_sdm(3),ind_dm_zone(6),ind_dm_omp(6),topo_s(3),
     &  ind_dm_socket(6), socket_topology(3), param_int(0:*), lok(*)

      REAL_E rop(*), rop_n(*), rop_n1(*), drodt(*)

      REAL_E param_real(0:*)

C Var loc 
      INTEGER_E nitrun,l,j,k,lij
#include "FastS/HPC_LAYER/LOC_VAR_DECLARATION.for"

      REAL_E c1,c2,c3, dtinv

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"

      !coeficient pour calcul gradient ordre2 !
      dtinv = 1./param_real(DTC)

      c1 =  1.5*dtinv
      c2 = -2.*dtinv
      c3 =  0.5*dtinv

      nitrun = -2

#include "FastS/HPC_LAYER/SIZE_MIN.for"
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
