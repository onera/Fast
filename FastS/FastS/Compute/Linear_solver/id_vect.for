c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine id_vect(param_int , ind_loop,
     &                   drodmd, krylov_out,  krylov_in)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    selection preconditionneur a droite pour Krylov
c
c     VAL
c
c     I/O
c_/    vectin    : vecteur kryloc IN
c_/    vectout   : vecteur kryloc OUT
c
c***********************************************************************     
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*)

      REAL_E     drodmd( param_int(NDIMDX), param_int(NEQ))
      REAL_E krylov_out( param_int(NDIMDX), param_int(NEQ)) 
      REAL_E  krylov_in( param_int(NDIMDX), param_int(NEQ))

C Var loc
      INTEGER_E k, j, lij, l,ltij,lt,lvo

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
             krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(l, 1)
             krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(l, 2)
             krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(l, 3)
             krylov_out(l, 4) = 0.
             krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(l, 5)
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
             krylov_out(l, 1) = (- drodmd(l, 1) + krylov_in(l, 1))
             krylov_out(l, 2) = (- drodmd(l, 2) + krylov_in(l, 2))
             krylov_out(l, 3) = (- drodmd(l, 3) + krylov_in(l, 3))
             krylov_out(l, 4) = 0.
             krylov_out(l, 5) = (- drodmd(l, 5) + krylov_in(l, 5))
             krylov_out(l, 6) = (- drodmd(l, 6) + krylov_in(l, 6))
#include   "FastS/Compute/loop_end.for"
        endif

      !!
      !!zone 3d
      !!
      ELSE

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
               krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(l, 1)
               krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(l, 2)
               krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(l, 3)
               krylov_out(l, 4) = - drodmd(l, 4) + krylov_in(l, 4)
               krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(l, 5)
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
               krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(l, 1)
               krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(l, 2)
               krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(l, 3)
               krylov_out(l, 4) = - drodmd(l, 4) + krylov_in(l, 4)
               krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(l, 5)
               krylov_out(l, 6) = - drodmd(l, 6) + krylov_in(l, 6)
#include   "FastS/Compute/loop_end.for"
        endif
      ENDIF

      end

