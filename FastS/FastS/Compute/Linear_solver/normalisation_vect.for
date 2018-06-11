      subroutine normalisation_vect(normL2, param_int, ind_loop, krylov)

      implicit none

#include "FastS/param_solver.h"
      
      INTEGER_E ind_loop(6), param_int(0:*)
      REAL_E normL2, krylov(param_int(NDIMDX), param_int(NEQ))

      INTEGER_E k, j, lij, l

#include "FastS/formule_param.h"

      do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)

            lij = inddm(ind_loop(1), j, k)

            do l = lij, lij + ind_loop(2) - ind_loop(1)

               krylov(l, 1) = krylov(l, 1) / normL2
               krylov(l, 2) = krylov(l, 2) / normL2
               krylov(l, 3) = krylov(l, 3) / normL2
               krylov(l, 4) = 0.0
               krylov(l, 5) = krylov(l, 5) / normL2
               
            enddo
         enddo
      enddo

      end
