      subroutine norm_vect( param_int, ind_loop, krylov, norm)

      implicit none

#include "FastS/param_solver.h"
      
      INTEGER_E ind_loop(6), param_int(0:*)
      REAL_E normL2, krylov(param_int(NDIMDX), param_int(NEQ)), norm

      INTEGER_E k, j, lij, l

#include "FastS/formule_param.h"

      !norm =0.
      do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3)-2, ind_loop(4)+2

            lij = inddm(ind_loop(1)-2, j, k)

            do l = lij, lij + ind_loop(2) - ind_loop(1) +4

               norm = norm + krylov(l, 1)**1
               norm = norm + krylov(l, 2)**1
               norm = norm + krylov(l, 3)**1
               norm = norm + krylov(l, 5)**1
               
            enddo
         enddo
      enddo

      end
