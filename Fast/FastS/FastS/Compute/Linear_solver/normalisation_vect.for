      subroutine normalisation_vect(normL2, param_int, ind_loop, krylov)

      implicit none

#include "FastS/param_solver.h"
      
      INTEGER_E ind_loop(6), param_int(0:*)
      REAL_E normL2, krylov(param_int(NDIMDX), param_int(NEQ))

C Var loc
      INTEGER_E k, j, lij, l,ltij,lt,lvo,i
      REAL_E norm_i

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      norm_i = 1./normL2

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) = krylov(l, 1) * norm_i
             krylov(l, 2) = krylov(l, 2) * norm_i
             krylov(l, 3) = krylov(l, 3) * norm_i
             krylov(l, 4) = 0.
             krylov(l, 5) = krylov(l, 5) * norm_i
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) = krylov(l, 1) * norm_i
             krylov(l, 2) = krylov(l, 2) * norm_i
             krylov(l, 3) = krylov(l, 3) * norm_i
             krylov(l, 4) = 0.
             krylov(l, 5) = krylov(l, 5) * norm_i
             krylov(l, 6) = krylov(l, 6) * norm_i
#include   "FastS/Compute/loop_end.for"
        endif

      !!
      !!zone 3d
      !!
      ELSE

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) = krylov(l, 1) * norm_i
             krylov(l, 2) = krylov(l, 2) * norm_i
             krylov(l, 3) = krylov(l, 3) * norm_i
             krylov(l, 4) = krylov(l, 4) * norm_i
             krylov(l, 5) = krylov(l, 5) * norm_i
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) = krylov(l, 1) * norm_i
             krylov(l, 2) = krylov(l, 2) * norm_i
             krylov(l, 3) = krylov(l, 3) * norm_i
             krylov(l, 4) = krylov(l, 4) * norm_i
             krylov(l, 5) = krylov(l, 5) * norm_i
             krylov(l, 6) = krylov(l, 6) * norm_i
#include   "FastS/Compute/loop_end.for"
        endif
      ENDIF

      end
