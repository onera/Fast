c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: Mehmet Demirci $
c***********************************************************************
      subroutine dpj_dpw_calc_meth(ndom, ithread,
     &                    param_int, param_real, param_int_eff,
     &                    ind_loop, xyz_ref,
     &                    rop, wig, xmut,
     &                    x, y, z,
     &                    dpCdp_dpW, dpClp_dpW,
     &                    ti, tj,tk, vol, venti, ventj, ventk,
     &                    cosAoA, sinAoA, surfinv)
c***********************************************************************
c_U   USER : Mehmet Demirci
c
c     ACT
c_A    Appel du calcul des effx explicites
c
c     VAL
c_V    gaz parfait monoespece
c_V    processeur domaine
c_V    steady/unsteady
c
c     INP
c_I    tijk     : vecteur normale aux facettes des mailles
c_I    ventijk     : vitesses d entrainement aux facettes preced.
c_I    qm,qp    : etats droit et gauche aux interfaces d une maille
c
c     LOC
c_L    eff      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    drodm    : increment de la solution
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ithread, ind_loop(6), param_int(0:*),
     & param_int_eff(0:*)


      REAL_E rop(*),xmut(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*),xyz_ref(*),
     & x(*),y(*),z(*),
     & dpCdp_dpW(*), dpClp_dpW(*),
     & cosAoA, sinAoA, surfinv

      REAL_E param_real(0:*)

C Var loc
      INTEGER_E kflu_loc

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

       !Si ALE et/ou SA, on force le passage par non ALE et NSLam
       !Si  Roe, on passe par AUSM pour avoir bon flux paroi

       kflu_loc = param_int(KFLUDOM)
       !if(param_int(KFLUDOM).eq.5) kflu_loc =1

       !write(*,'(a,7i7)')'optio=', option,ind_loop

      IF (kflu_loc.eq.1) THEN

          call dpJ_dpW_fluausm_select(ndom, ithread,
     &                        param_int, param_real,param_int_eff,
     &                        ind_loop, xyz_ref,
     &                        rop, wig,
     &                        x, y, z,
     &                        dpCdp_dpW, dpClp_dpW,
     &                        venti, ventj, ventk,
     &                        ti, tj, tk, vol, xmut,
     &                        cosAoA, sinAoA, surfinv)

      ELSEIF (kflu_loc.eq.2) THEN

c          call dpJ_dpW_flusenseur_init_select(ndom, ithread,
c     &                        param_int, param_real,param_int_eff,
c     &                        ind_loop, effort, xyz_ref,
c     &                        rop, flu, wig,
c     &                        x, y, z,
c     &                        venti, ventj, ventk,
c     &                        ti, tj, tk, vol, xmut)

      ELSEIF (kflu_loc.eq.5) THEN

c          call dpJ_dpW_fluroe_select(ndom, ithread,
c     &                        param_int, param_real,param_int_eff,
c     &                        ind_loop, effort, xyz_ref,
c     &                        rop, flu, wig,
c     &                        x, y, z,
c     &                        venti, ventj, ventk,
c     &                        ti, tj, tk, vol, xmut)
!!
!!!!
!!!!  ROE minmod: on ne passe plus par là, car flux paroi not equal flux roe
!!!!
!!
      ELSE

            write(*,*)'effort non code'
           call error('bceffort$',70,1)

      ENDIF
 
      end

