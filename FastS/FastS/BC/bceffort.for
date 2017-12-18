c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bceffort(ndom, ithread,
     &                    param_int, param_real, param_int_eff,
     &                    ind_loop, effort, xyz_ref,
     &                    rop, flu, wig, xmut,
     &                    x, y, z,
     &                    ti, tj,tk, vol, venti, ventj, ventk)
c***********************************************************************
c_U   USER : PECHIER
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


      REAL_E rop(*),xmut(*),flu(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*),effort(*),xyz_ref(*),
     & x(*),y(*),z(*)

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

          call eff_fluausm_select(ndom, ithread,
     &                        param_int, param_real,param_int_eff,
     &                        ind_loop, effort, xyz_ref,
     &                        rop, flu, wig,
     &                        x, y, z,
     &                        venti, ventj, ventk,
     &                        ti, tj, tk, vol, xmut)

      ELSEIF (kflu_loc.eq.2) THEN

          call eff_flusenseur_init_select(ndom, ithread,
     &                        param_int, param_real,param_int_eff,
     &                        ind_loop, effort, xyz_ref,
     &                        rop, flu, wig,
     &                        x, y, z,
     &                        venti, ventj, ventk,
     &                        ti, tj, tk, vol, xmut)

      ELSEIF (kflu_loc.eq.5) THEN

          call eff_fluroe_select(ndom, ithread,
     &                        param_int, param_real,param_int_eff,
     &                        ind_loop, effort, xyz_ref,
     &                        rop, flu, wig,
     &                        x, y, z,
     &                        venti, ventj, ventk,
     &                        ti, tj, tk, vol, xmut)

c         ndf    = param_int_eff(EFF_NDF)
c
c         pt_bc  = param_int(BC_NBBC + 1 + ndf)
c         bc_type= param_int(pt_bc + BC_TYPE)
c
c#include "FastS/BC/CL_correction_flu.for"

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

