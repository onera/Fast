      subroutine pre_bc(param_int, signe, ind_loop, krylov, rop)

      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*)
      REAL_E rop( param_int(NDIMDX), param_int(NEQ) )
      REAL_E krylov( param_int(NDIMDX), param_int(NEQ) )
      REAL_E signe

      INTEGER_E i,k, j, lij, l, lt,ltij,lvo

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return 

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) =  krylov(l, 1) + rop(l, 1)*signe
             krylov(l, 2) =  krylov(l, 2) + rop(l, 2)*signe
             krylov(l, 3) =  krylov(l, 3) + rop(l, 3)*signe
             krylov(l, 5) =  krylov(l, 5) + rop(l, 5)*signe
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) =  krylov(l, 1) + rop(l, 1)*signe
             krylov(l, 2) =  krylov(l, 2) + rop(l, 2)*signe
             krylov(l, 3) =  krylov(l, 3) + rop(l, 3)*signe
             krylov(l, 5) =  krylov(l, 5) + rop(l, 5)*signe
             krylov(l, 6) =  krylov(l, 6) + rop(l, 6)*signe
#include   "FastS/Compute/loop_end.for"
        endif

      !!
      !!zone 3d
      !!
      ELSE

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) =  krylov(l, 1) + rop(l, 1)*signe
             krylov(l, 2) =  krylov(l, 2) + rop(l, 2)*signe
             krylov(l, 3) =  krylov(l, 3) + rop(l, 3)*signe
             krylov(l, 4) =  krylov(l, 4) + rop(l, 4)*signe
             krylov(l, 5) =  krylov(l, 5) + rop(l, 5)*signe
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
             krylov(l, 1) =  krylov(l, 1) + rop(l, 1)*signe
             krylov(l, 2) =  krylov(l, 2) + rop(l, 2)*signe
             krylov(l, 3) =  krylov(l, 3) + rop(l, 3)*signe
             krylov(l, 4) =  krylov(l, 4) + rop(l, 4)*signe
             krylov(l, 5) =  krylov(l, 5) + rop(l, 5)*signe
             krylov(l, 6) =  krylov(l, 6) + rop(l, 6)*signe
#include   "FastS/Compute/loop_end.for"
        endif
      ENDIF

      end
