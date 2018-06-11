      subroutine bvbs_wall_inviscid_jacob(param_int, ind_loop,
     &                                    idir, neq_mtr,
     &                                    tijk, vect)

      implicit none

#include "FastS/param_solver.h"

      INTEGER_E param_int(0:*), ind_loop(6), neq_mtr, idir
      REAL_E tijk(param_int(NDIMDX_MTR), neq_mtr), vect(*)

      INTEGER_E k, j, lij, ltij, l, lt, v1, v2, v3, v4, v5, lp3,
     &     inc1, inc2, cellr, lijn,iref,i,ldjr
      REAL_E ncx, ncy, nsurf, dpfic1_dp1, dpfic2_dp2, dpfic2_dp3, 
     &     dpfic3_dp2, dpfic3_dp3, dpfic5_dp5
      

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      v1 = 0
      v2 = param_int(NDIMDX)
      v3 = 2 * param_int(NDIMDX)
      v4 = 3 * param_int(NDIMDX)
      v5 = 4 * param_int(NDIMDX)

      !write(*,*)'loop bc', ind_loop, idir
!     Utiliser shape_tab_mtr pour gerer 2d/3d


      IF (idir.eq.1) THEN

       iref = 2*ind_loop(2) + 1

          do 100 k = ind_loop(5), ind_loop(6)
          do 100 j = ind_loop(3), ind_loop(4)
          do 100 i = ind_loop(1), ind_loop(2)

            l    = inddm(  i             , j,  k ) 
            ldjr = inddm(  iref - i      , j,  k )
            lt   = indmtr( ind_loop(2)+1 , j,  k )

            ncx = tijk(lt, 1)
            ncy = tijk(lt, 2)
            nsurf = max(sqrt(ncx**2 + ncy**2), 1e-30)
            nsurf = 1. / nsurf
            ncx = ncx * nsurf
            ncy = ncy * nsurf

            dpfic1_dp1 = 1.
            dpfic2_dp2 = 1. - 2. * ncx**2
            dpfic2_dp3 = - 2. * ncx * ncy
            dpfic3_dp2 = - 2. * ncx * ncy
            dpfic3_dp3 = 1. - 2. * ncy**2
            dpfic5_dp5 = 1.

            !Premiere range de fictive
            vect(l + v1) = dpfic1_dp1 * vect(ldjr + v1 )
            vect(l + v2) = dpfic2_dp2 * vect(ldjr + v2 ) +
     &                     dpfic2_dp3 * vect(ldjr + v3 )
            vect(l + v3) = dpfic3_dp2 * vect(ldjr + v2 ) +
     &                     dpfic3_dp3 * vect(ldjr + v3 )
            vect(l + v4) = 0.
            vect(l + v5) = dpfic5_dp5 * vect(ldjr + v5 )
100    continue

      ELSEIF (idir.eq.2) THEN

       iref = 2*ind_loop(1) - 1

          do 120 k = ind_loop(5), ind_loop(6)
          do 120 j = ind_loop(3), ind_loop(4)
          do 120 i = ind_loop(1), ind_loop(2)

              l    = inddm(  i           , j, k ) 
              ldjr = inddm(  iref - i    , j, k )
              lt   = indmtr( ind_loop(1) , j, k )

            ncx = tijk(lt, 1)
            ncy = tijk(lt, 2)
            nsurf = max(sqrt(ncx**2 + ncy**2), 1e-30)
            nsurf = 1. / nsurf
            ncx = ncx * nsurf
            ncy = ncy * nsurf

            dpfic1_dp1 = 1.
            dpfic2_dp2 = 1. - 2. * ncx**2
            dpfic2_dp3 = - 2. * ncx * ncy
            dpfic3_dp2 = - 2. * ncx * ncy
            dpfic3_dp3 = 1. - 2. * ncy**2
            dpfic5_dp5 = 1.

            !Premiere range de fictive
            vect(l + v1) = dpfic1_dp1 * vect(ldjr + v1 )
            vect(l + v2) = dpfic2_dp2 * vect(ldjr + v2 ) +
     &                     dpfic2_dp3 * vect(ldjr + v3 )
            vect(l + v3) = dpfic3_dp2 * vect(ldjr + v2 ) +
     &                     dpfic3_dp3 * vect(ldjr + v3 )
            vect(l + v4) = 0.
            vect(l + v5) = dpfic5_dp5 * vect(ldjr + v5 )

120    continue

      else

        if (idir == 3) then
         j     = ind_loop(4)
         inc1  = param_int(NIJK)
         inc2  = 2 * param_int(NIJK)
         cellr = 1
        elseif (idir == 4) then
         j     = ind_loop(3)
         inc1  = - param_int(NIJK)
         inc2  = - 2 * param_int(NIJK)
         cellr = - 1
        endif

        do k = ind_loop(5), ind_loop(6)
         
         lij  = inddm(ind_loop(1), j, k)
         lijn = inddm(ind_loop(1), j + cellr, k)
         ltij = lijn - indmtr(ind_loop(1), j + cellr, k)
         
         do l = lij, lij + ind_loop(2) - ind_loop(1)

            lt = l - lij + lijn - ltij
            ncx = tijk(lt, 1)
            ncy = tijk(lt, 2)
            nsurf = max(sqrt(ncx**2 + ncy**2), 1e-30)
            nsurf = 1. / nsurf
            ncx = ncx * nsurf
            ncy = ncy * nsurf

            dpfic1_dp1 = 1.
            dpfic2_dp2 = 1. - 2. * ncx**2
            dpfic2_dp3 = - 2. * ncx * ncy
            dpfic3_dp2 = - 2. * ncx * ncy
            dpfic3_dp3 = 1. - 2. * ncy**2
            dpfic5_dp5 = 1.

            !Premiere range de fictive
            vect(l + v1) = dpfic1_dp1 * vect(l + v1 + inc1)
            vect(l + v2) = dpfic2_dp2 * vect(l + v2 + inc1) +
     &           dpfic2_dp3 * vect(l + v3 + inc1)
            vect(l + v3) = dpfic3_dp2 * vect(l + v2 + inc1) +
     &           dpfic3_dp3 * vect(l + v3 + inc1)
            vect(l + v4) = 0.
            vect(l + v5) = dpfic5_dp5 * vect(l + v5 + inc1)

           !Deuxieme range de fictive
            vect(l + v1 - inc1) = dpfic1_dp1 * vect(l + v1 + inc2)
            vect(l + v2 - inc1) = dpfic2_dp2 * vect(l + v2 + inc2) +
     &           dpfic2_dp3 * vect(l + v3 + inc2)
            vect(l + v3 - inc1) = dpfic3_dp2 * vect(l + v2 + inc2) +
     &           dpfic3_dp3 * vect(l + v3 + inc2)
            vect(l + v4 - inc1) = 0.
            vect(l + v5 - inc1) = dpfic5_dp5 * vect(l + v5 + inc2)

         enddo
      enddo
      endif

      end
