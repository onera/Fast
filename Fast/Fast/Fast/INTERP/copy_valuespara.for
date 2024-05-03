c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010)
c     $ $Revision: 56 $
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copy_valuespara(param_int, ind_loop,ind_loop_,rop,
     &                           stock,ind,sol,taille)
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

      INTEGER_E ind_loop(6), ind_loop_(6), param_int(0:*),ind,taille
      INTEGER_E sol

      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      !REAL_E stock(taille,param_int(NEQ))
      REAL_E stock(taille)
      !REAL_E stock(param_int(NDIMDX),param_int(NEQ))

C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
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
      !print*, "sol = ", sol

      nistk  =(ind_loop_(2)- ind_loop_(1))+1
      nistk2 =(ind_loop_(4)- ind_loop_(3))+1
      nistk3 =(ind_loop_(6)- ind_loop_(5))+1

      !print*, ind_loop(1)," ",ind_loop(2)
      !print*, ind_loop(3)," ",ind_loop(4)
      !print*, ind_loop(5)," ",ind_loop(6)

      !print*, "copy_valuespara"

      if (ind.eq.1) then  ! Stockage


         !do  ne=1,neq
         do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
               do  i = ind_loop(1), ind_loop(2)
                  do ne = 1, neq
                     
                     l  = inddm(i,j,k)
                     
                     cur  = (i - ind_loop_(1) + 1)
     &                    + (j - ind_loop_(3))*nistk
     &                    + (k - ind_loop_(5))*nistk*nistk2
                     lstk = neq*timelevels*(cur-1)
     &                       + neq*sol + ne
c                     lstk = neq*timelevels*cur
c     &                       + neq*sol + ne
c                     print*, lstk

                     stock(lstk) = rop(l,ne)
                     !if (ne.eq.1) then 
                     !print*, l, " ", rop(l,ne)
                     !endif
                     !rop(l,ne) = 100.
                     !if (ne.eq.1) then
                     !print*, rop(l,ne)
                     !endif
                     !print*, 'l= ', l," ",param_int(LEVEL)

                     !if(j==ind_loop(4).and.i==ind_loop(2) ) then
                     !print*,"stocke= ",rop(l,1),stock(lstk),l,lstk
                     !end if

                  enddo
                  !cur = cur+1
                  !print*, cur
               enddo
            enddo
         enddo
         !enddo


      elseif (ind.eq.2) then  ! Recuperation des valeurs

        !print*, "coucou recup valeurs"

        !do  ne=1,neq
        do  k = ind_loop(5), ind_loop(6)
           do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)
                 do ne = 1, neq

                    l  = inddm(i,j,k)

                     cur  = (i - ind_loop_(1) + 1)
     &                    + (j - ind_loop_(3))*nistk
     &                    + (k - ind_loop_(5))*nistk*nistk2
                     lstk = neq*timelevels*(cur-1)
     &                       + neq*sol + ne
c                     lstk = neq*timelevels*cur
c     &                       + neq*sol + ne

                    rop(l,ne) =  stock(lstk)

                    !if(j==ind_loop(4).and.i==ind_loop(2) ) then
                    !print*,"rop recup= ",rop(l,ne),k,lstk,ne
                    !end if

                 enddo
                 !cur = cur+1
                 !print*, cur
              enddo
           enddo
        enddo
c         enddo

c$$$       if (ind_loop(1)==-1 .or. ind_loop(1)==0) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = ind_loop(3), ind_loop(4)
c$$$                do  i = ind_loop(1), 0
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i+1,j,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$
c$$$       else if (ind_loop(2)==param_int(NIJK) .or.
c$$$     & ind_loop(2)==param_int(NIJK)-1) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = ind_loop(3), ind_loop(4)
c$$$                do  i = param_int(NIJK)-1 , ind_loop(2)
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i-1,j,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$
c$$$
c$$$       else if (ind_loop(3)==-1 .or. ind_loop(3)==0) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = ind_loop(3), 0
c$$$                do  i = ind_loop(1) , ind_loop(2)
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i,j+1,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$
c$$$
c$$$       else if (ind_loop(4)==param_int(NIJK+1) .or.
c$$$     & ind_loop(4)==param_int(NIJK+1)-1) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = param_int(NIJK+1)-1, ind_loop(4)
c$$$                do  i = ind_loop(1) , ind_loop(2)
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i,j-1,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$
c$$$
c$$$        end if




      end if







      end
