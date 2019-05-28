c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine bvbs_periodique(idir,lrhs, param_int, ind_loop, rop)
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

      INTEGER_E idir,lrhs, ind_loop(6), param_int(0:*)

      REAL_E rop( param_int(NDIMDX), param_int(NEQ) )

C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,inc1,ne,lij
      REAL_E c1

#include "FastS/formule_param.h"

      ! mise a zero si tableau drodm du RHS implicit
      if(lrhs.eq.1) then
        c1= 0.
      else !extrapolation ordre 0 si variable conservative ou primitive
        c1= 1.
      endif

      if(idir.eq.1.or.idir.eq.2) then

        inc1     = param_int(IJKV)
        if(idir.eq.2) inc1 = -param_int(IJKV)

         do  ne=1,param_int(NEQ)
         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         
            do  i = ind_loop(1), ind_loop(2)

              ldjr     =  inddm(inc1 + i ,j,k)
              l         = inddm(i        ,j,k)

              rop(l,ne) = rop(ldjr,ne)*c1
            enddo
         enddo
         enddo
         enddo

      elseif(idir.eq.3.or.idir.eq.4) then

        inc1     = param_int(IJKV+1)*param_int(NIJK)
        if(idir.eq.4) inc1 = -inc1

         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         
            lij  = inddm(ind_loop(1) , j, k)

#ifdef _OPENMP4
!$OMP simd 
#else
CDIR$ IVDEP
#endif
            do  l = lij, lij + ind_loop(2) - ind_loop(1)

              ldjr      = l + inc1

              rop(l,1) = rop(ldjr,1)*c1
              rop(l,2) = rop(ldjr,2)*c1
              rop(l,3) = rop(ldjr,3)*c1
              rop(l,4) = rop(ldjr,4)*c1
              rop(l,5) = rop(ldjr,5)*c1
            enddo
           if(param_int(NEQ).ge.6) then
              do  ne=6,param_int(NEQ)
#ifdef _OPENMP4
!$OMP simd
#else
CDIR$ IVDEP
#endif
              do  l = lij, lij + ind_loop(2) - ind_loop(1)

                ldjr      = l + inc1
 
                rop(l,ne) = rop(ldjr,ne)*c1
              enddo
             enddo
           endif

         enddo
         enddo

      else

        inc1     = param_int(IJKV+2)*param_int(NIJK)*param_int(NIJK+1)
        if(idir.eq.6) inc1  = -inc1

         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)

            lij  = inddm(ind_loop(1) , j, k)

#ifdef _OPENMP4
!$OMP simd
#else
CDIR$ IVDEP
#endif
            do  l = lij, lij + ind_loop(2) - ind_loop(1)

              ldjr      = l + inc1

              rop(l,1) = rop(ldjr,1)*c1
              rop(l,2) = rop(ldjr,2)*c1
              rop(l,3) = rop(ldjr,3)*c1
              rop(l,4) = rop(ldjr,4)*c1
              rop(l,5) = rop(ldjr,5)*c1
            enddo

           if(param_int(NEQ).ge.6) then
             do  ne=6,param_int(NEQ)
#ifdef _OPENMP4
!$OMP simd
#else
CDIR$ IVDEP
#endif
              do  l = lij, lij + ind_loop(2) - ind_loop(1)

                ldjr      = l + inc1
 
                rop(l,ne) = rop(ldjr,ne)*c1
              enddo
             enddo
           endif
         enddo
         enddo

      endif

c      write(*,*)'idir=', idir,inc1,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop=', ind_loop

      end
