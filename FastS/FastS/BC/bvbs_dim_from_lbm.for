c***********************************************************************
      subroutine bvbs_dim_from_lbm(idir, lrhs, param_int, param_real,
     &                             ind_loop, rop)
c***********************************************************************
c
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"


      INTEGER_E lrhs, ind_loop(6), param_int(0:*)
      INTEGER_E idir
      REAL_E param_real(0:*)

      REAL_E rop( param_int(NDIMDX), param_int(NEQ) )
C Var local
      INTEGER_E l,i,j,k, l_prev,l_pre2
      INTEGER_E startk,endk,startj,endj,starti,endi
      REAL_E c1
      REAL_E pp, rgp, gam, tinf, roinf, pp1,pp2
      REAL_E scalingro, scalingu, scalingt, c0, dt, dx
      REAL_E psolver, puser, p0solver, p0user, rolbm
      REAL_E uinf, adv, tfluclbm, tdelta, tnsextrap

#include "FastS/formule_param.h"

      !Recuperation des grandeurs de reference pour l'adim
      !print*, 'bc dim ns'
      dt    = param_real( DTC )
      !print*, 'dt = ',dt
      gam   = param_real( GAMMA )
      !print*, 'gam = ',gam
      rgp   = param_real( CVINF )*(gam-1.)  !Cv(gama-1)= R (gas parfait)
      !print*, 'rgp = ',rgp
      tinf  = param_real( TINF )
      !print*, 'tinf = ',tinf
      roinf = param_real( ROINF )
      !print*, 'roinf = ',roinf
      p0user  = param_real( PINF )
      !print*, 'p0user = ',p0user
      uinf    = param_real( VINF )
      !print*, 'uinf = ',uinf
      !scalingro = 1.
      scalingro = roinf
      c0 = sqrt(gam*rgp*tinf)
      scalingu = sqrt(3.)*c0
      !scalingt = 300.101375
      scalingt = tinf
      !print*, tinf
      !print*, roinf
      dx    = dt*c0*sqrt(3.)

      ! Mise a zero si tableau drodm du RHS implicit
C      if(lhrs.eq.1) then
C        c1= 0.
C      else !extrapolation ordre 0 si variable conservative ou primitive
C        c1=1.
C      endif

c      startk = min(1-param_int(NIJK+3),ind_loop(5))
c      endk   = max(2+param_int(IJKV+2),ind_loop(6))
c      startj = min(1-param_int(NIJK+3),ind_loop(3))
c      endj   = max(2+param_int(IJKV+1),ind_loop(4))

      starti = ind_loop(1)
      endi   = ind_loop(2)
      startj = ind_loop(3)
      endj   = ind_loop(4)
      startk = ind_loop(5)
      endk   = ind_loop(6)

      ! CURRENT VERSION
      ! if (ind_loop(1)==1 - param_int(NIJK+3)) then
      !     starti = ind_loop(1)
      !     endi = ind_loop(2)
      !     if (ind_loop(5)==1) then
      !        startk = 1 - param_int(NIJK+3)
      !     else
      !        startk = ind_loop(5)
      !     endif
      !     if (ind_loop(6)==param_int(IJKV+2)) then
      !        endk = param_int(IJKV+2) + param_int(NIJK+3)
      !     else
      !        endk = ind_loop(6)
      !     endif
      !     if (ind_loop(3)==1) then
      !        startj = 1 - param_int(NIJK+3)
      !     else
      !        startj = ind_loop(3)
      !     endif
      !     if (ind_loop(4)==param_int(IJKV+1)) then
      !        endj = param_int(IJKV+1) + param_int(NIJK+3)
      !     else
      !        endj = ind_loop(4)
      !     endif
      ! elseif (ind_loop(2)==param_int(IJKV) + param_int(NIJK+3)) then
      !     starti = ind_loop(1)
      !     endi = ind_loop(2)
      !     if (ind_loop(5)==1) then
      !        startk = 1 - param_int(NIJK+3)
      !     else
      !        startk = ind_loop(5)
      !     endif
      !     if (ind_loop(6)==param_int(IJKV+2)) then
      !        endk = param_int(IJKV+2) + param_int(NIJK+3)
      !     else
      !        endk = ind_loop(6)
      !     endif
      !     if (ind_loop(3)==1) then
      !        startj = 1 - param_int(NIJK+3)
      !     else
      !        startj = ind_loop(3)
      !     endif
      !     if (ind_loop(4)==param_int(IJKV+1)) then
      !        endj = param_int(IJKV+1) + param_int(NIJK+3)
      !     else
      !        endj = ind_loop(4)
      !     endif
      ! else
      !   !startk = ind_loop(5)
      !   !endk = ind_loop(6)
      !   if (ind_loop(5)==1) then
      !      startk = 1 - param_int(NIJK+3)
      !   else
      !      startk = ind_loop(5)
      !   endif
      !   if (ind_loop(6)==param_int(IJKV+2)) then
      !      endk = param_int(IJKV+2) + param_int(NIJK+3)
      !   else
      !      endk = ind_loop(6)
      !   endif
      !   startj = ind_loop(3)
      !   endj = ind_loop(4)
      !   starti = ind_loop(1)
      !   endi = ind_loop(2)
      ! endif

      ! if (ind_loop(5)==1) then
      !   startk = 1 - param_int(NIJK+3)
      ! else
      !   startk = ind_loop(5)
      ! endif
      ! if (ind_loop(6)==param_int(IJKV+2)) then
      !   endk = param_int(IJKV+2) + param_int(NIJK+3)
      ! else
      !   endk = ind_loop(6)
      ! endif
      ! if (ind_loop(3)==1) then
      !   startj = 1 - param_int(NIJK+3)
      ! else
      !   startj = ind_loop(3)
      ! endif
      ! if (ind_loop(4)==param_int(IJKV+1)) then
      !   endj = param_int(IJKV+1) + param_int(NIJK+3)
      ! else
      !   endj = ind_loop(4)
      ! endif
      ! if (ind_loop(1)==1) then
      !   starti = 1 - param_int(NIJK+3)
      ! else
      !   starti = ind_loop(1)
      ! endif
      ! if (ind_loop(2)==param_int(IJKV)) then
      !   endi = param_int(IJKV) + param_int(NIJK+3)
      ! else
      !   endi = ind_loop(2)
      ! endif

      ! print*, 'Boundary NS'
      ! print*, 'imin', starti, 'imax', endi
      ! print*, 'jmin', startj, 'jmax', endj
      ! print*, 'kmin', startk, 'kmax', endk
      ! print*, ''

      !do k = ind_loop(5), ind_loop(6)
      do k = startk, endk
        do j = startj, endj
        !do j = ind_loop(3), ind_loop(4)
          do i = starti, endi
          !do i = ind_loop(1), ind_loop(2)

             l = inddm(i,j,k)
             l_prev = inddm(i-1,j,k)
             l_pre2 = inddm(i-2,j,k)

             ! On reconstruit T a partir de la pression
             rolbm = rop(l,1)/roinf
             pp = roinf*gam*rgp*tinf*( -0.4/1.4 + rolbm )
             tfluclbm = pp/(rop(l,1)*rgp)
             ! rop(l,5) = pp/(rop(l,1)*rgp)

             ! Extrapolation
             ! rop(l,5) = rop(l_prev,5)
             ! rop(l,5) = rop(l_pre2,5)+2*(rop(l_prev,5)-rop(l_pre2,5))

             ! Advection
             ! rop(l,5) = rop(l_prev,5)
             ! adv = -uinf*dt/dx*(rop(l,5)-rop(l_prev,5))
             ! rop(l,5) = rop(l,5)+adv

             !Enrichissement de la temperature LBM
             tnsextrap = rop(l_prev,5)
             tdelta = tnsextrap-tfluclbm!-tinf
             !rop(l,5) = tinf + tdelta

              ! if (idir==1.or.idir==3.or.idir==4) then
              !   rop(l,5) = pp/(rop(l,1)*rgp)
              ! endif


          end do
        end do
      end do

      end subroutine
