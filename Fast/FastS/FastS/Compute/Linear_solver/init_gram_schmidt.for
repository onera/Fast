      subroutine init_gram_schmidt(param_int, ind_loop, drodm,
     &     krylov, normL2)

      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*)
      REAL_E drodm(param_int(NDIMDX), param_int(NEQ)),
     &     krylov(param_int(NDIMDX), param_int(NEQ)), normL2

      INTEGER_E k, j, lij, l
      REAL_E tmp1, tmp2, tmp3, tmp5

#include "FastS/formule_param.h"

      if (ind_loop(1) .GT. ind_loop(2)) return
      if (ind_loop(3) .GT. ind_loop(4)) return
      if (ind_loop(5) .GT. ind_loop(6)) return

      normL2 = 0.

!     V0
      do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)

            lij = inddm(ind_loop(1), j, k)

            do l = lij, lij + ind_loop(2) - ind_loop(1)

               tmp1 = drodm(l, 1)
               tmp2 = drodm(l, 2)
               tmp3 = drodm(l, 3)
               tmp5 = drodm(l, 5)

               krylov(l, 1) = tmp1
               krylov(l, 2) = tmp2
               krylov(l, 3) = tmp3
               krylov(l, 4) = 0.0
               krylov(l, 5) = tmp5

               normL2 = normL2 + tmp1 * tmp1 + tmp2 * tmp2 +
     &              tmp3 * tmp3 + tmp5 * tmp5

            enddo
         enddo
      enddo

      end

