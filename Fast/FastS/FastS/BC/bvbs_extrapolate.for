c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine bvbs_extrapolate(idir,lrhs, eq_deb,
     &                            param_int, ind_loop, nut_inf, rop)
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

#include "FastS/param_solver.h"

      INTEGER_E idir,lrhs, eq_deb, ind_loop(6), param_int(0:*)

      REAL_E rop( param_int(NDIMDX), param_int(NEQ) ), nut_inf
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,inc1,ne,lij
      REAL_E c1, vmax

#include "FastS/formule_param.h"

      ! mise a zero si tableau drodm du RHS implicit
      if(lrhs.eq.1) then
        c1= 0.
      else !extrapolation ordre 0 si variable conservative ou primitive
        c1= 1.
      endif

      vmax = -1.e30
      if(eq_deb.eq.6) vmax = nut_inf

      if(idir.eq.1.or.idir.eq.2) then

        inc1     = 1
        if(idir.eq.2) inc1 = ind_loop(1)-1

         do  ne=eq_deb,param_int(NEQ)
         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         
         ldjr = inddm(inc1,j,k)
            do  i = ind_loop(1), ind_loop(2)

              l         = inddm(i,j,k)

              rop(l,ne) = max( vmax, rop(ldjr,ne) ) *c1
            enddo
         enddo
         enddo
         enddo

      elseif(idir.eq.3.or.idir.eq.4) then

        inc1     = 1
        if(idir.eq.4) inc1 = ind_loop(3)-1

         do  ne=eq_deb,param_int(NEQ)
         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         
            do  i = ind_loop(1), ind_loop(2)

              ldjr      = inddm(i, inc1, k)
              l         = inddm(i,    j, k)

              rop(l,ne) = max( vmax, rop(ldjr,ne) ) *c1
            enddo
         enddo
         enddo
         enddo

      else

       inc1     = 1
        if(idir.eq.6) inc1 = ind_loop(5)-1

         do  ne=eq_deb,param_int(NEQ)
         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         
            do  i = ind_loop(1), ind_loop(2)

              ldjr      = inddm(i, j, inc1)
              l         = inddm(i, j, k   )

              rop(l,ne) = max( vmax, rop(ldjr,ne) ) *c1
            enddo
         enddo
         enddo
         enddo

      endif

c      write(*,*)'idir=', idir,inc1,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop=', ind_loop

      end
