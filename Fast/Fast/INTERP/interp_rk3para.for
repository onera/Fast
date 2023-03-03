c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010)
c     $ $Revision: 56 $
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine interp_rk3para( param_int, ind_loop, ind_loop_,
     &                           rop, stock, nstep, nitrun, taille )
c***********************************************************************
c_U   USER : PECHIER ********SUBROUTINE FONCTIONNANT SEULE**************
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

      INTEGER_E ind_loop(6),ind_loop_(6),param_int(0:*),taille,nstep
      INTEGER_E nitrun

      REAL_E stock(taille)
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))

C Var local

      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk,ind
      INTEGER_E nistk2, nistk3
      INTEGER_E lstk_n, lstk_np1, lstk_nm1, timelevels, order
      INTEGER_E cur, sol_0, sol_1, sol_2
      REAL_E coeff1, coeff2, coeff3
      REAL_E ro,u,v,w,t

#include "Fast/formule_param.h"

      order = 1
      timelevels = 3
      neq = 5

      cur = 0

      sol_0 = 0
      sol_1 = 1
      sol_2 = 2

      nistk  =(ind_loop_(2)- ind_loop_(1))+1
      nistk2 =(ind_loop_(4)- ind_loop_(3))+1
      nistk3 =(ind_loop_(6)- ind_loop_(5))+1

      !print*, "nstep = ", nstep
      !print*, "nitrun = ", nitrun, order
      if (order==0) then
        coeff1 = 0.
        coeff2 = 1.
        coeff3 = 0.
      else if (order==1) then
        if (nstep==1.and.nitrun==1) then
           coeff1 = 0.5
           coeff2 = 0.5
           coeff3 = 0.
        else if (nstep==1.and.nitrun.ge.2) then
           coeff1 = 0.5!1.
           coeff2 = 0.5!0.25
           coeff3 = 0.!-0.25
        else if (nstep==2.and.nitrun==1) then
           coeff1 = 0.2113248654051877
           coeff2 = 0.7886751345948123
           coeff3 = 0.
        else if (nstep==2.and.nitrun.ge.2) then
           coeff1 = 0.2113248654051877!1.
           coeff2 = 0.7886751345948123!0.394337567
           coeff3 = 0.!-0.394337567
        end if
      else if (order==2) then
        if (nstep==1.and.nitrun==1) then
           coeff1 = 0.!0.5!0.
           coeff2 = 1.!0.5!1.
           coeff3 = 0.
        else if (nstep==1.and.nitrun.ge.2) then
           coeff1 = 0.!1.25  !1.25
           coeff2 = 1.!0.125 !-0.375
           coeff3 = 0.!-0.375 !0.125
        else if (nstep==2.and.nitrun==1) then
           coeff1 = 0.2113248654051877
           coeff2 = 0.7886751345948123
           coeff3 = 0.
        else if (nstep==2.and.nitrun.ge.2) then
           coeff1 = 1.166666666666667
           coeff2 = 0.3110042339640727
           coeff3 = -0.47767090063073964
        end if
      end if

      if (order.le.1) then

       do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
           do  i = ind_loop(1), ind_loop(2)

           l  = inddm(i,j,k)

           cur  = (i - ind_loop_(1) + 1)
     &          + (j - ind_loop_(3))*nistk
     &          + (k - ind_loop_(5))*nistk*nistk2

           lstk_n   = neq*timelevels*(cur-1) + neq*sol_0
           lstk_np1 = neq*timelevels*(cur-1) + neq*sol_1 

           rop(l,1) = coeff1*stock(lstk_n +1) +coeff2*stock(lstk_np1 +1)
           rop(l,2) = coeff1*stock(lstk_n +2) +coeff2*stock(lstk_np1 +2)
           rop(l,3) = coeff1*stock(lstk_n +3) +coeff2*stock(lstk_np1 +3)
           rop(l,4) = coeff1*stock(lstk_n +4) +coeff2*stock(lstk_np1 +4)
           rop(l,5) = coeff1*stock(lstk_n +5) +coeff2*stock(lstk_np1 +5)

           end do
         end do
       end do

      else 
      !print*, "Je passe dans la fonction d'interp"
      !print*, "coef1 = ", coeff1, " coef2 = ",coeff2," coef3 = ",coeff3

      do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
          do  i = ind_loop(1), ind_loop(2)

               l  = inddm(i,j,k)

               cur  = (i - ind_loop_(1) + 1)
     &              + (j - ind_loop_(3))*nistk
     &              + (k - ind_loop_(5))*nistk*nistk2

               lstk_n   = neq*timelevels*(cur-1)
     &                    + neq*sol_0 + 1
               lstk_np1 = neq*timelevels*(cur-1)
     &                    + neq*sol_1 + 1
               lstk_nm1 = neq*timelevels*(cur-1)
     &                    + neq*sol_2 + 1

c               lstk_n   = neq*timelevels*cur
c     &                    + neq*sol_0 + 1
c               lstk_np1 = neq*timelevels*cur
c     &                    + neq*sol_1 + 1
c               lstk_nm1 = neq*timelevels*cur
c     &                    + neq*sol_2 + 1

c               print*, "lstk_n = ", lstk_n, " lstk_np1 = ", lstk_np1,
c     &                 " lstk_nm1 = ", lstk_nm1

               ro = coeff1*stock(lstk_n)+coeff2*stock(lstk_np1)
     &              + coeff3*stock(lstk_nm1)

      write(*,*)'ro',stock(lstk_n),stock(lstk_np1),stock(lstk_nm1),k,j,i
               !if(j==ind_loop(4).and.i==ind_loop(2)) then
               !print*,l,lstk_n,lstk_np1,lstk_nm1
               !print*, stock(lstk_n),stock(lstk_np1),stock(lstk_nm1)
               !print*,"ro apres interp = ", ro
               !print*,"rop interp =", rop(l,3)
               !print*,"rop interp =", rop(l,4)
               !print*,"rop interp =", rop(l,5)
               !endif
               lstk_n   = neq*timelevels*(cur-1)
     &                    + neq*sol_0 + 2
               lstk_np1 = neq*timelevels*(cur-1)
     &                    + neq*sol_1 + 2
               lstk_nm1 = neq*timelevels*(cur-1)
     &                    + neq*sol_2 + 2

c               lstk_n   = neq*timelevels*cur
c     &                    + neq*sol_0 + 2
c               lstk_np1 = neq*timelevels*cur
c     &                    + neq*sol_1 + 2
c               lstk_nm1 = neq*timelevels*cur
c     &                    + neq*sol_2 + 2

      write(*,*)'u ',stock(lstk_n),stock(lstk_np1),stock(lstk_nm1)
               u = coeff1*stock(lstk_n)+coeff2*stock(lstk_np1)
     &             + coeff3*stock(lstk_nm1)

               lstk_n   = neq*timelevels*(cur-1)
     &                    + neq*sol_0 + 3
               lstk_np1 = neq*timelevels*(cur-1)
     &                    + neq*sol_1 + 3
               lstk_nm1 = neq*timelevels*(cur-1)
     &                    + neq*sol_2 + 3

c               lstk_n   = neq*timelevels*cur
c     &                    + neq*sol_0 + 3
c               lstk_np1 = neq*timelevels*cur
c     &                    + neq*sol_1 + 3
c               lstk_nm1 = neq*timelevels*cur
c     &                    + neq*sol_2 + 3

      write(*,*)'v ',stock(lstk_n),stock(lstk_np1),stock(lstk_nm1)
               v = coeff1*stock(lstk_n)+coeff2*stock(lstk_np1)
     &             + coeff3*stock(lstk_nm1)

               lstk_n   = neq*timelevels*(cur-1)
     &                    + neq*sol_0 + 4
               lstk_np1 = neq*timelevels*(cur-1)
     &                    + neq*sol_1 + 4
               lstk_nm1 = neq*timelevels*(cur-1)
     &                    + neq*sol_2 + 4

c               lstk_n   = neq*timelevels*cur
c     &                    + neq*sol_0 + 4
c               lstk_np1 = neq*timelevels*cur
c     &                    + neq*sol_1 + 4
c               lstk_nm1 = neq*timelevels*cur
c     &                    + neq*sol_2 + 4

      write(*,*)'w ',stock(lstk_n),stock(lstk_np1),stock(lstk_nm1)
               w = coeff1*stock(lstk_n)+coeff2*stock(lstk_np1)
     &             + coeff3*stock(lstk_nm1)

               lstk_n   = neq*timelevels*(cur-1)
     &                    + neq*sol_0 + 5
               lstk_np1 = neq*timelevels*(cur-1)
     &                    + neq*sol_1 + 5
               lstk_nm1 = neq*timelevels*(cur-1)
     &                    + neq*sol_2 + 5

c               lstk_n   = neq*timelevels*cur
c     &                    + neq*sol_0 + 5
c               lstk_np1 = neq*timelevels*cur
c     &                    + neq*sol_1 + 5
c               lstk_nm1 = neq*timelevels*cur
c     &                    + neq*sol_2 + 5

      write(*,*)'T ',stock(lstk_n),stock(lstk_np1),stock(lstk_nm1)

               t = coeff1*stock(lstk_n)+coeff2*stock(lstk_np1)
     &             + coeff3*stock(lstk_nm1)

               !if(j==ind_loop(4).and.i==ind_loop(2)) then
               !print*,"rop interp =", rop(l,1),l,lstk_np1
               !print*,"rop interp =", rop(l,2)
               !print*,"rop interp =", rop(l,3)
               !print*,"rop interp =", rop(l,4)
               !print*,"rop interp =", rop(l,5)
               !endif
               rop(l,1) = ro
               rop(l,2) = u
               rop(l,3) = v
               rop(l,4) = w
               rop(l,5) = t

               !if(j==ind_loop(4).and.i==ind_loop(2)) then
               !print*,"rop interp =", rop(l,1)
               !print*,"rop interp =", rop(l,2)
               !print*,"rop interp =", rop(l,3)
               !print*,"rop interp =", rop(l,4)
               !print*,"rop interp =", rop(l,5)
               !endif
               !cur = cur + 1

              end do
            end do
         end do

         endif


         end
