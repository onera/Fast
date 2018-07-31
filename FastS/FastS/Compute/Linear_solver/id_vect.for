c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine id_vect(param_int , ind_loop,
     &                   drodmd, krylov_out,  krylov_in, size)
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

      INTEGER_E ind_loop(6), param_int(0:*),size

      REAL_E     drodmd( param_int(NDIMDX), param_int(NEQ))
      REAL_E krylov_out( param_int(NDIMDX), param_int(NEQ)) 
      REAL_E  krylov_in( size             , param_int(NEQ))

C Var loc
      INTEGER_E i, k, j, l, i_2,j_2,k_2,v1,v2,v3,v4,ls,
     & v5,i_size,j_size,indlu

#include "FastS/formule_param.h"

      indlu(i_2,j_2,k_2) = 1 + (i_2+v1) + (j_2+v2)*v4 + (k_2+v3)*v5

      if (param_int(NB_RELAX) .GE. 2) then

         i_size = ind_loop(2) - ind_loop(1) + 1 +
     &        2 * param_int(NIJK + 3) !taille de la fenetre + ghostcells
         j_size = ind_loop(4) - ind_loop(3) + 1 +
     &        2 * param_int(NIJK + 3)

         v1 = param_int(NIJK+3) - ind_loop(1)
         v2 = param_int(NIJK+3) - ind_loop(3)
         v3 = param_int(NIJK+4) - ind_loop(5)
         v4 = i_size
         v5 = i_size * j_size

      else

         v1 = param_int(NIJK+3) - 1
         v2 = param_int(NIJK+3) - 1
         v3 = param_int(NIJK+4) - 1
         v4 = param_int(NIJK)
         v5 = param_int(NIJK) * param_int(NIJK+1)

      endif

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, 1)
              ls = indlu(i, j, 1)

              krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(ls, 1)
              krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(ls, 2)
              krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(ls, 3)
              krylov_out(l, 4) = 0.
              krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(ls, 5)

           enddo
           enddo
        else
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, 1)
              ls = indlu(i, j, 1)

              krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(ls, 1)
              krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(ls, 2)
              krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(ls, 3)
              krylov_out(l, 4) = 0.
              krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(ls, 5)
              krylov_out(l, 6) = - drodmd(l, 6) + krylov_in(ls, 6)

           enddo
           enddo
        endif

      !!
      !!zone 3d
      !!
      ELSE

        if(param_int(NEQ).eq.5) then
           do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, k)
              ls = indlu(i, j, k)

              krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(ls, 1)
              krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(ls, 2)
              krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(ls, 3)
              krylov_out(l, 4) = - drodmd(l, 4) + krylov_in(ls, 4)
              krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(ls, 5)

           enddo
           enddo
           enddo
        else
           do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, k)
              ls = indlu(i, j, k)

              krylov_out(l, 1) = - drodmd(l, 1) + krylov_in(ls, 1)
              krylov_out(l, 2) = - drodmd(l, 2) + krylov_in(ls, 2)
              krylov_out(l, 3) = - drodmd(l, 3) + krylov_in(ls, 3)
              krylov_out(l, 4) = - drodmd(l, 4) + krylov_in(ls, 4)
              krylov_out(l, 5) = - drodmd(l, 5) + krylov_in(ls, 5)
              krylov_out(l, 6) = - drodmd(l, 6) + krylov_in(ls, 6)

           enddo
           enddo
           enddo
        endif
      ENDIF

      end

