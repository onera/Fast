c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine prod_mat_vect(param_int, ind_loop, krylov, vecty,
     &     drodm, num_krylov)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    selection preconditionneur a droite pour Krylov
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

      INTEGER_E param_int(0:*), ind_loop(6), num_krylov

      REAL_E drodm(  param_int(NDIMDX) * param_int(NEQ) )
      REAL_E krylov( param_int(NDIMDX) * param_int(NEQ), num_krylov)
      REAL_E vecty(num_krylov - 1)

!var loc
      INTEGER_E k,j,lij, l, col, v1, v2, v3,v4,v5,v6,ltij,lt,lvo

      REAL_E sum_prod1, sum_prod2, sum_prod3, sum_prod5

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v6 = 5*param_int(NDIMDX)


      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

       if(param_int(NEQ).eq.5) then

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

          ! mise a zero pour calcul somme
#include  "FastS/Compute/loopI_begin.for"
            drodm(l + v1) = 0.
            drodm(l + v2) = 0.
            drodm(l + v3) = 0.
            drodm(l + v5) = 0.
          enddo

           do col = 1, num_krylov - 1

#include    "FastS/Compute/loopI_begin.for"
             drodm(l+ v1)= drodm(l+ v1) + krylov(l+ v1, col)*vecty(col)
             drodm(l+ v2)= drodm(l+ v2) + krylov(l+ v2, col)*vecty(col)
             drodm(l+ v3)= drodm(l+ v3) + krylov(l+ v3, col)*vecty(col)
             drodm(l+ v5)= drodm(l+ v5) + krylov(l+ v5, col)*vecty(col)
            enddo
           enddo
c
        enddo
        enddo

       else

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

          ! mise a zero pour calcul somme
#include  "FastS/Compute/loopI_begin.for"
            drodm(l + v1) = 0.
            drodm(l + v2) = 0.
            drodm(l + v3) = 0.
            drodm(l + v5) = 0.
            drodm(l + v6) = 0.
          enddo

           do col = 1, num_krylov - 1

#include    "FastS/Compute/loopI_begin.for"
             drodm(l+ v1)= drodm(l+ v1) + krylov(l+ v1, col)*vecty(col)
             drodm(l+ v2)= drodm(l+ v2) + krylov(l+ v2, col)*vecty(col)
             drodm(l+ v3)= drodm(l+ v3) + krylov(l+ v3, col)*vecty(col)
             drodm(l+ v5)= drodm(l+ v5) + krylov(l+ v5, col)*vecty(col)
             drodm(l+ v6)= drodm(l+ v6) + krylov(l+ v6, col)*vecty(col)
            enddo
           enddo
c
        enddo
        enddo
       endif ! 5-6 eq

      !3d
      ELSE  

       if(param_int(NEQ).eq.5) then

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

          ! mise a zero pour calcul somme
#include  "FastS/Compute/loopI_begin.for"
            drodm(l + v1) = 0.
            drodm(l + v2) = 0.
            drodm(l + v3) = 0.
            drodm(l + v4) = 0.
            drodm(l + v5) = 0.
          enddo

           do col = 1, num_krylov - 1

#include    "FastS/Compute/loopI_begin.for"
             drodm(l+ v1)= drodm(l+ v1) + krylov(l+ v1, col)*vecty(col)
             drodm(l+ v2)= drodm(l+ v2) + krylov(l+ v2, col)*vecty(col)
             drodm(l+ v3)= drodm(l+ v3) + krylov(l+ v3, col)*vecty(col)
             drodm(l+ v4)= drodm(l+ v4) + krylov(l+ v4, col)*vecty(col)
             drodm(l+ v5)= drodm(l+ v5) + krylov(l+ v5, col)*vecty(col)
            enddo
           enddo
c
        enddo
        enddo

       else

        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)

          ! mise a zero pour calcul somme
#include  "FastS/Compute/loopI_begin.for"
            drodm(l + v1) = 0.
            drodm(l + v2) = 0.
            drodm(l + v3) = 0.
            drodm(l + v4) = 0.
            drodm(l + v5) = 0.
            drodm(l + v6) = 0.
          enddo

           do col = 1, num_krylov - 1

#include    "FastS/Compute/loopI_begin.for"
             drodm(l+ v1)= drodm(l+ v1) + krylov(l+ v1, col)*vecty(col)
             drodm(l+ v2)= drodm(l+ v2) + krylov(l+ v2, col)*vecty(col)
             drodm(l+ v3)= drodm(l+ v3) + krylov(l+ v3, col)*vecty(col)
             drodm(l+ v4)= drodm(l+ v4) + krylov(l+ v4, col)*vecty(col)
             drodm(l+ v5)= drodm(l+ v5) + krylov(l+ v5, col)*vecty(col)
             drodm(l+ v6)= drodm(l+ v6) + krylov(l+ v6, col)*vecty(col)
            enddo
           enddo
c
        enddo
        enddo
       endif ! 5-6 eq

      ENDIF

      end
