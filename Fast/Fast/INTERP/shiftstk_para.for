c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010)
c     $ $Revision: 56 $
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine shiftstk_para(param_int, ind_loop,ind_loop_,
     &                           stock,sol_i, sol_d, taille)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    extrapolation ordere zero cell fictive
c
c     VAL
c_V    Optimisation NEC
c
c     COM
c***********************************************************************
      implicit none

#include "Fast/param_solver.h"

      INTEGER_E ind_loop(6), ind_loop_(6), param_int(0:*),taille
      INTEGER_E sol_i, sol_d

      REAL_E stock(taille)

C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk
      INTEGER_E lstk_i, lstk_d
      INTEGER_E nistk2,nistk3,cycl
      INTEGER_E timelevels, cur
      REAL_E c1


#include "Fast/formule_param.h"

      neq=5
      timelevels = 3
      cur = 0




      !print*,"param_int(NEQ)= ",param_int(NEQ)
      !print*, "ind= ",ind
      !print*, pos
      !print*, "taille stock = ", taille
      !print*, "sol initial = ", sol_i
      !print*, "sol destination = ", sol_d

      nistk  =(ind_loop_(2)- ind_loop_(1))+1
      nistk2 =(ind_loop_(4)- ind_loop_(3))+1
      nistk3 =(ind_loop_(6)- ind_loop_(5))+1

      !print*, ind_loop(1)," ",ind_loop(2)
      !print*, ind_loop(3)," ",ind_loop(4)
      !print*, ind_loop(5)," ",ind_loop(6)

      !print*, "coucou shift"

      !do  ne=1,neq
      do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)
               do ne = 1, neq

                  l  = inddm(i,j,k)

                  cur  = (i - ind_loop_(1) + 1)
     &                 + (j - ind_loop_(3))*nistk
     &                 + (k - ind_loop_(5))*nistk*nistk2

                  lstk_i = neq*timelevels*(cur-1)
     &                     + neq*sol_i + ne
c                  lstk_i = neq*timelevels*cur
c     &                     + neq*sol_i + ne
                  !print*, lstk_i

                  lstk_d = neq*timelevels*(cur-1)
     &                     + neq*sol_d + ne
c                  lstk_d = neq*timelevels*cur
c     &                     + neq*sol_d + ne
                  !print*, lstk_d

                  stock(lstk_d) = stock(lstk_i)
                  !if (ne.eq.1) then
                  !print*, l, " ", rop(l,ne)
                  !endif
                  !rop(l,ne) = 100.
                  !if (ne.eq.1) then
                  !print*, rop(l,ne)
                  !endif
                  !print*, 'l= ', l," ",param_int(LEVEL)

                  !if(j==ind_loop(4).and.i==ind_loop(2) ) then
                  !print*,"mv stk = ",stock(lstk_d),stock(lstk_i),cur,ne
                  !end if

               enddo

               !cur = cur+1
               !print*, cur

            enddo
         enddo
      enddo
       !enddo

      end
