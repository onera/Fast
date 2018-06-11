c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine LDURelaxScal(param_int, param_real, ind_loop,
     &                        step, nobl, norm,
     &                        vectin, vectout, rop, coe, ti, tj, ssor)
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

      INTEGER_E ind_loop(6), param_int(0:*), nb_relax, nobl

      REAL_E param_real(0:*), ti(*), tj(*)
      REAL_E    ssor( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E     rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E  vectin( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E vectout( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E     coe( param_int(NDIMDX) , param_int(NEQ_COE) )

      REAL_E  norm

C Var loc
      INTEGER_E indcell(nobl), indinm(nobl), indjnm(nobl),
     &     indsmtr(nobl), indin1(nobl), indjn1(nobl)

      INTEGER_E step, inci, incj, incimj, incimjmtr, im, jm, v1mtr,
     &     v2mtr, nsommax, nsommin, nsom, nptobl,k,lij,ltij,l,lt,lvo,
     &     nptoblinm, nptobljnm, ibeg, iend, indclbeg, indmtbeg,
     &     i, j, np, ncl, nmtr, nmtrpi, nmtrpj, npjnm, nclpj,
     &     npinm, nclpi, nptoblin1, npjn1, nclmj, npin1, nclmi, 
     &     nptobljn1

      REAL_E  nx, ny, nz, surf, uu, vv, ww, ee, hh, g1ee, g1uu, g1vv,
     &     g1ww, gam, gam1, ce, coef, ext1, ext2, ext3, ext4, ext5, cp,
     &     vsn, ttt, dw1, dw2, dw3, dw4, dw5, diag, rspecA, rspecB, c,
     &     valuescal,norm_1

      REAL_E am11, am12, am13, am14, am15, am21, am22, am23, am24, am25,
     &     am31, am32, am33, am34, am35, am41, am42, am43, am44, am45, 
     &     am51, am52, am53, am54, am55

      REAL_E bm11, bm12, bm13, bm14, bm15, bm21, bm22, bm23, bm24, bm25,
     &     bm31, bm32, bm33, bm34, bm35, bm41, bm42, bm43, bm44, bm45, 
     &     bm51, bm52, bm53, bm54, bm55

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      gam  = param_real(GAMMA)
      ce   = 1.5
      gam1 = gam - 1.
      cp   = gam * param_real(CVINF)

      inci      = 1
      incj      = param_int(NIJK)
      incimj    = inci - incj
      incimjmtr = param_int(NIJK_MTR) - param_int(NIJK_MTR + 1)

      im        = ind_loop(2) + 1
      jm        = ind_loop(4) + 1

      v1mtr = 0
      v2mtr = param_int(NDIMDX_MTR)

      nsommax = im + jm - 2
      nsommin = 2


      norm_1 = 1./norm

      if (step == 1) then
#include   "FastS/Compute/loop_begin.for"
            vectout(l,1) = vectin(l,1)*norm_1
            vectout(l,2) = vectin(l,2)*norm_1
            vectout(l,3) = vectin(l,3)*norm_1
            vectout(l,4) = vectin(l,4)*norm_1
            vectout(l,5) = vectin(l,5)*norm_1
               ssor(l,1) = 0.
               ssor(l,2) = 0.
               ssor(l,3) = 0.
               ssor(l,4) = 0.
               ssor(l,5) = 0.
#include   "FastS/Compute/loop_end.for"
       else
#include   "FastS/Compute/loop_begin.for"
            vectout(l,1) = vectin(l,1) + ssor(l,1)
            vectout(l,2) = vectin(l,2) + ssor(l,2)
            vectout(l,3) = vectin(l,3) + ssor(l,3)
            vectout(l,4) = vectin(l,4) + ssor(l,4)
            vectout(l,5) = vectin(l,5) + ssor(l,5)
               ssor(l,1) = 0.
               ssor(l,2) = 0.
               ssor(l,3) = 0.
               ssor(l,4) = 0.
               ssor(l,5) = 0.
#include   "FastS/Compute/loop_end.for"
      endif
      


      if (mod(step,2) == 1) then
!     faire LD
         do nsom = nsommin, nsommax

            nptobl = 0

            nptoblinm  = 0
            nptobljnm  = 0

            ibeg = max(1, nsom - jm + 1)
            iend = min(im - 1, nsom - 1)
            indclbeg = inddm(ibeg, nsom - ibeg, 1)
            indmtbeg = indmtr(ibeg, nsom - ibeg, 1)

            do i = ibeg, iend

               j = nsom - i
               nptobl = nptobl + 1
               indcell(nptobl) = indclbeg + (i - ibeg) * incimj
               indsmtr(nptobl) = indmtbeg + (i - ibeg) *
     &              incimjmtr

               if (i .NE. (im - 1)) then

                  nptoblinm = nptoblinm + 1
                  indinm(nptoblinm) = nptobl

               endif

               if (j .NE. (jm - 1)) then

                  nptobljnm = nptobljnm + 1
                  indjnm(nptobljnm) = nptobl

               endif

            enddo               !fin boucle i

            do np = 1, nptobl

               ncl = indcell(np)
               nmtr = indsmtr(np)
               nmtrpi = nmtr + param_int(NIJK_MTR)
               nmtrpj = nmtr + param_int(NIJK_MTR + 1)

               uu = rop(ncl, 2)
               vv = rop(ncl, 3)
               ee = 0.5 * (uu**2 + vv**2)
               hh = cp * rop(ncl, 5) + ee
               c = sqrt(gam1 * (hh - ee))

               nx = 0.5 * (tj(nmtr + v1mtr) + tj(nmtrpj + v1mtr))
               ny = 0.5 * (tj(nmtr + v2mtr) + tj(nmtrpj + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               vsn = uu * nx + vv * ny
               rspecB = abs(vsn) + c * surf

               nx = 0.5 * (ti(nmtr + v1mtr) + ti(nmtrpi + v1mtr))
               ny = 0.5 * (ti(nmtr + v2mtr) + ti(nmtrpi + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               vsn = uu * nx + vv * ny
               rspecA = abs(vsn) + c * surf

               diag = ce * coe(ncl, 1) * (rspecA + rspecB) + 1.
               diag = 1. / diag

               vectout(ncl, :) = vectout(ncl, :) * diag

            enddo               !fin boucle np

            do npjnm = 1, nptobljnm

               np = indjnm(npjnm)
               ncl = indcell(np)
               nclpj = ncl + incj
               nmtr = indsmtr(np)
               nmtrpj = indsmtr(np) + param_int(NIJK_MTR + 1)

               uu = rop(ncl, 2)
               vv = rop(ncl, 3)
               ww = 0.
               ee = 0.5 * (uu**2 + vv**2)
               hh = cp * rop(ncl, 5) + ee
               c = sqrt(gam1 * (hh - ee))

               g1ee = gam1 * ee
               g1uu = gam1 * uu
               g1vv = gam1 * vv
               g1ww = gam1 * ww

               nx = 0.5 * (tj(nmtr + v1mtr) + tj(nmtrpj + v1mtr))
               ny = 0.5 * (tj(nmtr + v2mtr) + tj(nmtrpj + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               nz = 0.
               vsn = uu * nx + vv * ny
               rspecB = abs(vsn) + c * surf

               bm11 = rspecB
               bm12 = nx
               bm13 = ny
               bm14 = nz
               bm15 = 0.

               bm21 = g1ee * nx - uu * vsn
               bm22 = (2. - gam) * uu * nx + vsn + rspecB
               bm23 = - g1vv * nx + uu * ny
               bm24 = - g1ww * nx + uu * nz
               bm25 = gam1 * nx

               bm31 = g1ee * ny - vv * vsn
               bm32 = - g1uu * ny + vv * nx
               bm33 = (2. - gam) * vv * ny + vsn + rspecB
               bm34 = - g1ww * ny+vv * nz
               bm35 = gam1 * ny

               bm41 = g1ee * nz - ww * vsn
               bm42 = - g1uu * nz + ww * nx
               bm43 = - g1vv * nz + ww * ny
               bm44 = (2. - gam) * ww * nz + vsn + rspecB
               bm45 = gam1 * nz

               bm51 = (g1ee - hh) * vsn
               bm52 = hh * nx - g1uu * vsn
               bm53 = hh * ny - g1vv * vsn
               bm54 = hh * nz - g1ww * vsn
               bm55 = gam * vsn + rspecB

               ttt = ce * 0.5 * coe(nclpj, 1)

               dw1 = ttt*(bm11 * vectout(ncl, 1) + bm12 * vectout(ncl,2)
     &              + bm13 * vectout(ncl, 3) + bm14 * vectout(ncl, 4)
     &              + bm15 * vectout(ncl, 5))

               dw2 = ttt*(bm21 * vectout(ncl, 1) + bm22 * vectout(ncl,2)
     &              + bm23 * vectout(ncl, 3) + bm24 * vectout(ncl, 4)
     &              + bm25 * vectout(ncl, 5))

               dw3 = ttt*(bm31 * vectout(ncl, 1) + bm32 * vectout(ncl,2)
     &              + bm33 * vectout(ncl, 3) + bm34 * vectout(ncl, 4)
     &              + bm35 * vectout(ncl, 5))

               dw4 = ttt*(bm41 * vectout(ncl, 1) + bm42 * vectout(ncl,2)
     &              + bm43 * vectout(ncl, 3) + bm44 * vectout(ncl, 4)
     &              + bm45 * vectout(ncl, 5))

               dw5 = ttt*(bm51 * vectout(ncl, 1) + bm52 * vectout(ncl,2)
     &              + bm53 * vectout(ncl, 3) + bm54 * vectout(ncl, 4)
     &              + bm55 * vectout(ncl, 5))

               vectout(nclpj, 1) = vectout(nclpj, 1) + dw1
               vectout(nclpj, 2) = vectout(nclpj, 2) + dw2
               vectout(nclpj, 3) = vectout(nclpj, 3) + dw3
               vectout(nclpj, 4) = vectout(nclpj, 4) + dw4
               vectout(nclpj, 5) = vectout(nclpj, 5) + dw5

               ssor(nclpj, 1) = ssor(nclpj, 1) + dw1
               ssor(nclpj, 2) = ssor(nclpj, 2) + dw2
               ssor(nclpj, 3) = ssor(nclpj, 3) + dw3
               ssor(nclpj, 4) = ssor(nclpj, 4) + dw4
               ssor(nclpj, 5) = ssor(nclpj, 5) + dw5

            enddo               !fin boucle npjnm

            do npinm = 1, nptoblinm

               np = indinm(npinm)
               ncl = indcell(np)
               nclpi = ncl + inci
               nmtr = indsmtr(np)
               nmtrpi = nmtr + param_int(NIJK_MTR)

               uu = rop(ncl, 2)
               vv = rop(ncl, 3)
               ww = 0.
               ee = 0.5 * (uu**2 + vv**2)
               hh = cp * rop(ncl, 5) + ee
               c = sqrt(gam1 * (hh - ee))

               g1ee = gam1 * ee
               g1uu = gam1 * uu
               g1vv = gam1 * vv
               g1ww = gam1 * ww

               nx = 0.5 * (ti(nmtr + v1mtr) + ti(nmtrpi + v1mtr))
               ny = 0.5 * (ti(nmtr + v2mtr) + ti(nmtrpi + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               nz = 0.
               vsn = uu * nx + vv * ny
               rspecA = abs(vsn) + c * surf

               am11 = rspecA
               am12 = nx
               am13 = ny
               am14 = nz
               am15 = 0.

               am21 = g1ee * nx - uu * vsn
               am22 = (2. - gam) * uu * nx + vsn + rspecA
               am23 = - g1vv * nx + uu * ny
               am24 = - g1ww * nx + uu * nz
               am25 = gam1 * nx

               am31 = g1ee * ny - vv * vsn
               am32 = - g1uu * ny + vv * nx
               am33 = (2. - gam) * vv * ny + vsn + rspecA
               am34 = - g1ww * ny + vv * nz
               am35 = gam1 * ny

               am41 = g1ee * nz - ww * vsn
               am42 = - g1uu * nz + ww * nx
               am43 = - g1vv * nz + ww * ny
               am44 = (2. - gam) * ww * nz + vsn + rspecA
               am45 = gam1 * nz

               am51 = (g1ee - hh) * vsn
               am52 = hh * nx - g1uu * vsn
               am53 = hh * ny - g1vv * vsn
               am54 = hh * nz - g1ww * vsn
               am55 = gam * vsn + rspecA

               ttt = ce * 0.5 * coe(nclpi, 1)

               dw1 = ttt*(am11 * vectout(ncl, 1) + am12 * vectout(ncl,2)
     &              + am13 * vectout(ncl, 3) + am14 * vectout(ncl, 4)
     &              + am15 * vectout(ncl, 5))

               dw2 = ttt*(am21 * vectout(ncl, 1) + am22 * vectout(ncl,2)
     &              + am23 * vectout(ncl, 3) + am24 * vectout(ncl, 4)
     &              + am25 * vectout(ncl, 5))

               dw3 = ttt*(am31 * vectout(ncl, 1) + am32 * vectout(ncl,2)
     &              + am33 * vectout(ncl, 3) + am34 * vectout(ncl, 4)
     &              + am35 * vectout(ncl, 5))

               dw4 = ttt*(am41 * vectout(ncl, 1) + am42 * vectout(ncl,2)
     &              + am43 * vectout(ncl, 3) + am44 * vectout(ncl, 4)
     &              + am45 * vectout(ncl, 5))

               dw5 = ttt*(am51 * vectout(ncl, 1) + am52 * vectout(ncl,2)
     &              + am53 * vectout(ncl, 3) + am54 * vectout(ncl, 4)
     &              + am55 * vectout(ncl, 5))

               vectout(nclpi, 1) = vectout(nclpi, 1) + dw1
               vectout(nclpi, 2) = vectout(nclpi, 2) + dw2
               vectout(nclpi, 3) = vectout(nclpi, 3) + dw3
               vectout(nclpi, 4) = vectout(nclpi, 4) + dw4
               vectout(nclpi, 5) = vectout(nclpi, 5) + dw5

               ssor(nclpi, 1) = ssor(nclpi, 1) + dw1
               ssor(nclpi, 2) = ssor(nclpi, 2) + dw2
               ssor(nclpi, 3) = ssor(nclpi, 3) + dw3
               ssor(nclpi, 4) = ssor(nclpi, 4) + dw4
               ssor(nclpi, 5) = ssor(nclpi, 5) + dw5

            enddo               !fin boucle npinm

         enddo                  !fin boucle nsom

      else

         do nsom = nsommax, nsommin, - 1

            nptobl = 0
            nptoblin1 = 0
            nptobljn1 = 0

            ibeg = max(1, nsom - jm + 1)
            iend = min(im - 1, nsom - 1)
            indclbeg = inddm(ibeg, nsom - ibeg, 1)
            indmtbeg = indmtr(ibeg, nsom - ibeg, 1)

            do i = ibeg, iend

               j = nsom - i
               nptobl = nptobl + 1
               indcell(nptobl) = indclbeg + (i - ibeg) * incimj
               indsmtr(nptobl) = indmtbeg + (i - ibeg) *
     &              incimjmtr

               if (i .NE. 1) then

                  nptoblin1 = nptoblin1 + 1
                  indin1(nptoblin1) = nptobl

               endif

               if (j .NE. 1) then

                  nptobljn1 = nptobljn1 + 1
                  indjn1(nptobljn1) = nptobl

               endif

            enddo               !fin boucle i

            do np = 1, nptobl

               ncl = indcell(np)
               nmtr = indsmtr(np)
               nmtrpi = indsmtr(np) + param_int(NIJK_MTR)
               nmtrpj = indsmtr(np) + param_int(NIJK_MTR + 1)

               uu = rop(ncl, 2)
               vv = rop(ncl, 3)
               ee = 0.5 * (uu**2 + vv**2)
               hh = cp * rop(ncl, 5) + ee
               c = sqrt(gam1 * (hh - ee))

               nx = 0.5 * (tj(nmtr + v1mtr) + tj(nmtrpj + v1mtr))
               ny = 0.5 * (tj(nmtr + v2mtr) + tj(nmtrpj + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               vsn = uu * nx + vv * ny
               rspecB = abs(vsn) + c * surf

               nx = 0.5 * (ti(nmtr + v1mtr) + ti(nmtrpi + v1mtr))
               ny = 0.5 * (ti(nmtr + v2mtr) + ti(nmtrpi + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               vsn = uu * nx + vv * ny
               rspecA = abs(vsn) + c * surf

               diag = ce * coe(ncl, 1) * (rspecA + rspecB) + 1.
               diag = 1. / diag

               vectout(ncl, :) = vectout(ncl, :) * diag

            enddo               !fin boucle np

            do npjn1 = 1, nptobljn1

               np = indjn1(npjn1)
               ncl = indcell(np)
               nclmj = ncl - incj
               nmtr = indsmtr(np)
               nmtrpj = nmtr + param_int(NIJK_MTR + 1)

               uu = rop(ncl, 2)
               vv = rop(ncl, 3)
               ww = 0.
               ee = 0.5 * (uu**2 + vv**2)
               hh = cp * rop(ncl, 5) + ee
               c = sqrt(gam1 * (hh - ee))

               g1ee = gam1 * ee
               g1uu = gam1 * uu
               g1vv = gam1 * vv
               g1ww = gam1 * ww

               nx = 0.5 * (tj(nmtr + v1mtr) + tj(nmtrpj + v1mtr))
               ny = 0.5 * (tj(nmtr + v2mtr) + tj(nmtrpj + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               nz = 0.
               vsn = uu * nx + vv * ny
               rspecB = abs(vsn) + c * surf

               bm11 = 0. - rspecB
               bm12 = nx
               bm13 = ny
               bm14 = nz
               bm15 = 0.

               bm21 = g1ee * nx - uu * vsn
               bm22 = (2. - gam) * uu * nx + vsn - rspecB
               bm23 = - g1vv * nx + uu * ny
               bm24 = - g1ww * nx + uu * nz
               bm25 = gam1 * nx

               bm31 = g1ee * ny - vv * vsn
               bm32 = - g1uu * ny + vv * nx
               bm33 = (2. - gam) * vv * ny + vsn  - rspecB
               bm34 = - g1ww * ny + vv * nz
               bm35 = gam1 * ny

               bm41 = g1ee * nz - ww * vsn
               bm42 = - g1uu * nz + ww * nx
               bm43 = - g1vv * nz + ww * ny
               bm44 = (2. - gam) * ww * nz + vsn - rspecB
               bm45 = gam1 * nz

               bm51 = (g1ee - hh) * vsn
               bm52 = hh * nx - g1uu * vsn
               bm53 = hh * ny - g1vv * vsn
               bm54 = hh * nz - g1ww * vsn
               bm55 = gam * vsn - rspecB

               ttt = ce * 0.5 * coe(nclmj, 1)

               dw1 = - ttt*(bm11 * vectout(ncl,1) + bm12 *vectout(ncl,2)
     &              + bm13 * vectout(ncl, 3) + bm14 * vectout(ncl, 4)
     &              + bm15 * vectout(ncl, 5))

               dw2 = - ttt*(bm21 * vectout(ncl,1) + bm22 *vectout(ncl,2)
     &              + bm23 * vectout(ncl, 3) + bm24 * vectout(ncl, 4)
     &              + bm25 * vectout(ncl, 5))

               dw3 = - ttt*(bm31 * vectout(ncl,1) + bm32 *vectout(ncl,2)
     &              + bm33 * vectout(ncl, 3) + bm34 * vectout(ncl, 4)
     &              + bm35 * vectout(ncl, 5))

               dw4 = - ttt*(bm41 * vectout(ncl,1) + bm42 *vectout(ncl,2)
     &              + bm43 * vectout(ncl, 3) + bm44 * vectout(ncl, 4)
     &              + bm45 * vectout(ncl, 5))

               dw5 = - ttt*(bm51 * vectout(ncl,1) + bm52 *vectout(ncl,2)
     &              + bm53 * vectout(ncl, 3) + bm54 * vectout(ncl, 4)
     &              + bm55 * vectout(ncl, 5))

               vectout(nclmj, 1) = vectout(nclmj, 1) + dw1
               vectout(nclmj, 2) = vectout(nclmj, 2) + dw2
               vectout(nclmj, 3) = vectout(nclmj, 3) + dw3
               vectout(nclmj, 4) = vectout(nclmj, 4) + dw4
               vectout(nclmj, 5) = vectout(nclmj, 5) + dw5

               ssor(nclmj, 1) = ssor(nclmj, 1) + dw1
               ssor(nclmj, 2) = ssor(nclmj, 2) + dw2
               ssor(nclmj, 3) = ssor(nclmj, 3) + dw3
               ssor(nclmj, 4) = ssor(nclmj, 4) + dw4
               ssor(nclmj, 5) = ssor(nclmj, 5) + dw5

            enddo               !fin boucle npjn1

            do npin1 = 1, nptoblin1

               np = indin1(npin1)
               ncl = indcell(np)
               nclmi = ncl - inci
               nmtr = indsmtr(np)
               nmtrpi = nmtr + param_int(NIJK_MTR)

               uu = rop(ncl, 2)
               vv = rop(ncl, 3)
               ww = 0.
               ee = 0.5 * (uu**2 + vv**2)
               hh = cp * rop(ncl, 5) + ee
               c = sqrt(gam1 * (hh - ee))

               g1ee = gam1 * ee
               g1uu = gam1 * uu
               g1vv = gam1 * vv
               g1ww = gam1 * ww

               nx = 0.5 * (ti(nmtr + v1mtr) + ti(nmtrpi + v1mtr))
               ny = 0.5 * (ti(nmtr + v2mtr) + ti(nmtrpi + v2mtr))
               surf = sqrt(nx**2 + ny**2)
               nz = 0.
               vsn = uu * nx + vv * ny
               rspecA = abs(vsn) + c * surf

               am11 = 0. - rspecA
               am12 = nx
               am13 = ny
               am14 = nz
               am15 = 0.

               am21 = g1ee * nx - uu * vsn
               am22 = (2. - gam) * uu * nx + vsn - rspecA
               am23 = - g1vv * nx + uu * ny
               am24 = - g1ww * nx + uu * nz
               am25 = gam1 * nx

               am31 = g1ee * ny - vv * vsn
               am32 = - g1uu * ny + vv * nx
               am33 = (2. - gam) * vv * ny + vsn - rspecA
               am34 = - g1ww * ny + vv * nz
               am35 = gam1 * ny

               am41 = g1ee * nz - ww * vsn
               am42 = - g1uu * nz + ww * nx
               am43 = - g1vv * nz + ww * ny
               am44 = (2. - gam) * ww * nz + vsn - rspecA
               am45 = gam1 * nz

               am51 = (g1ee - hh) * vsn
               am52 = hh * nx - g1uu * vsn
               am53 = hh * ny - g1vv * vsn
               am54 = hh * nz - g1ww * vsn
               am55 = gam * vsn - rspecA

               ttt = ce * 0.5 * coe(nclmi, 1)

               dw1 = - ttt*(am11 * vectout(ncl,1) + am12 *vectout(ncl,2)
     &              + am13 * vectout(ncl, 3) + am14 * vectout(ncl, 4)
     &              + am15 * vectout(ncl, 5))

               dw2 = - ttt*(am21 * vectout(ncl,1) + am22 *vectout(ncl,2)
     &              + am23 * vectout(ncl, 3) + am24 * vectout(ncl, 4)
     &              + am25 * vectout(ncl, 5))

               dw3 = - ttt*(am31 * vectout(ncl,1) + am32 *vectout(ncl,2)
     &              + am33 * vectout(ncl, 3) + am34 * vectout(ncl, 4)
     &              + am35 * vectout(ncl, 5))

               dw4 = - ttt*(am41 * vectout(ncl,1) + am42 *vectout(ncl,2)
     &              + am43 * vectout(ncl, 3) + am44 * vectout(ncl, 4)
     &              + am45 * vectout(ncl, 5))

               dw5 = - ttt*(am51 * vectout(ncl,1) + am52 *vectout(ncl,2)
     &              + am53 * vectout(ncl, 3) + am54 * vectout(ncl, 4)
     &              + am55 * vectout(ncl, 5))

               vectout(nclmi, 1) = vectout(nclmi, 1) + dw1
               vectout(nclmi, 2) = vectout(nclmi, 2) + dw2
               vectout(nclmi, 3) = vectout(nclmi, 3) + dw3
               vectout(nclmi, 4) = vectout(nclmi, 4) + dw4
               vectout(nclmi, 5) = vectout(nclmi, 5) + dw5

               ssor(nclmi, 1) = ssor(nclmi, 1) + dw1
               ssor(nclmi, 2) = ssor(nclmi, 2) + dw2
               ssor(nclmi, 3) = ssor(nclmi, 3) + dw3
               ssor(nclmi, 4) = ssor(nclmi, 4) + dw4
               ssor(nclmi, 5) = ssor(nclmi, 5) + dw5

            enddo               !fin boucle npin1

         enddo                  !fin boucle nsom

      endif

      end
