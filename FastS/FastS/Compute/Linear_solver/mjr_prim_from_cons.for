c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine mjr_prim_from_cons( param_int, param_real,
     &                               visco,sa_real,
     &                               ind_loop,
     &                               rop_1, rop, drodm_out)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    modifie variable primitive a partir increment conservatif
c
c     VAL
c
c     I/O
c_/    vect1,2  value :  IN
c_/    normL2         : OUT
c
c***********************************************************************²
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*)

      REAL_E visco(5),sa_real(3)
      REAL_E rop_1(param_int(NDIMDX), param_int(NEQ)), 
     &     rop(param_int(NDIMDX), param_int(NEQ)), 
     &     drodm_out(param_int(NDIMDX), param_int(NEQ)), param_real(0:*)

      INTEGER_E k, j, lij, l,v1, v2, v3,v4,v5,v6,ltij,lt,lvo

      REAL_E cvinv, cvinv2, ro_old, u_old, v_old, w_old,t_old, roe_old,
     &     r_1, anulam, nu_old, cmus1, temp01, coesut,ratiom

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      cvinv = 1. / param_real(CVINF)
      cvinv2 = 0.5 * cvinv

      cmus1  = visco(5)
      temp01 = 1./visco(4)
      coesut = visco(3) * (1.+cmus1*temp01)
      ratiom = sa_real(SA_RATIOM)

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

       if(param_int(NEQ).eq.5) then

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

#include  "FastS/Compute/loopI_begin.for"

#include  "FastS/Compute/LU/mjr_newton_2d.for"
          enddo

        enddo
        enddo

       else

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

#include  "FastS/Compute/loopI_begin.for"

#include  "FastS/Compute/LU/mjr_newton_2d_SA.for"
c           if (j.le.80.and.l-lij.eq.100)
c     &         write(*,*)ro_old*nu_old , drodm(l,6), rop_1(l,6)
          enddo
        enddo
        enddo
       endif ! 5-6 eq

      !3d
      ELSE  

       if(param_int(NEQ).eq.5) then

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

#include  "FastS/Compute/loopI_begin.for"

#include  "FastS/Compute/LU/mjr_newton.for"
          enddo

        enddo
        enddo

       else

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

#include  "FastS/Compute/loopI_begin.for"

#include  "FastS/Compute/LU/mjr_newton_SA.for"
          enddo
        enddo
        enddo
       endif ! 5-6 eq

      ENDIF

      end
