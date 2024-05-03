c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine bvbs_periodique_azimuthal(idir,lrhs, 
     &                                     param_int,ind_loop, 
     &                                     rop, data_per)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    periodicite geometrique (turbomachine)
c
c     VAL
c_V    Optimisation NEC
c
c     COM
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E idir,lrhs, ind_loop(6), param_int(0:*)

      REAL_E data_per(*)
      REAL_E rop( param_int(NDIMDX), param_int(NEQ) )

C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,inc1,ne,lij
      INTEGER_E dir_axe
      REAL_E c1
      REAL_E pi,angle,cosa,sina
      REAL_E ropr2,ropr3,ropr4
      REAL_E, DIMENSION(3,3) :: rot

#include "FastS/formule_param.h"

      ! matrice de rotation
      if    (data_per(2).eq.0. .and. data_per(3).eq.0.) then
         angle   = data_per(1)
         dir_axe = 1
      elseif(data_per(1).eq.0. .and. data_per(3).eq.0.) then
         angle   = data_per(2)
         dir_axe = 2
      elseif(data_per(1).eq.0. .and. data_per(2).eq.0.) then
         angle   = data_per(3)
         dir_axe = 3
      else
         WRITE(*,*) 'Warning: BC periodic: bad axis definition.'
         angle = 0.
         dir_axe = 1
      endif
     
      !correction angle
      angle =angle*(-1.)

      pi   = acos(-1.)
      cosa = cos(pi*angle/180.)
      sina = sin(pi*angle/180.)
      
      select case(dir_axe)
        
        case(1)
           
           rot(1,1) = 1.
           rot(1,2) = 0.
           rot(1,3) = 0.
           rot(2,1) = 0.
           rot(2,2) = cosa
           rot(2,3) = -sina
           rot(3,1) = 0.
           rot(3,2) = sina
           rot(3,3) = cosa
           
        case(2)
           
           rot(1,1) = cosa
           rot(1,2) = -sina
           rot(1,3) = 0.
           rot(2,1) = 0.
           rot(2,2) = 1.
           rot(2,3) = 0.
           rot(3,1) = sina
           rot(3,2) = cosa
           rot(3,3) = 0.  
           
        case(3)
           
           rot(1,1) = cosa
           rot(1,2) = -sina
           rot(1,3) = 0.
           rot(2,1) = sina
           rot(2,2) = cosa
           rot(2,3) = 0.
           rot(3,1) = 0.
           rot(3,2) = 0.
           rot(3,3) = 1.
        
        end select
        
      ! mise a zero si tableau drodm du RHS implicit
      if(lrhs.eq.1) then
        c1= 0.
      else !extrapolation ordre 0 si variable conservative ou primitive
        c1= 1.
      endif

      if(idir.eq.1.or.idir.eq.2) then

        do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)
        
              l         = inddm(i        ,j,k)
              
              ropr2 = rop(l,2)*rot(1,1)+rop(l,3)*rot(1,2)+
     &                rop(l,4)*rot(1,3)
              ropr3 = rop(l,2)*rot(2,1)+rop(l,3)*rot(2,2)+
     &                rop(l,4)*rot(2,3)
              ropr4 = rop(l,2)*rot(3,1)+rop(l,3)*rot(3,2)+
     &                rop(l,4)*rot(3,3)
            
              rop(l,1) = rop(l,1)*c1
              rop(l,2) = ropr2*c1
              rop(l,3) = ropr3*c1
              rop(l,4) = ropr4*c1
              rop(l,5) = rop(l,5)*c1
                
             enddo
          enddo
        enddo

      elseif(idir.eq.3.or.idir.eq.4) then

         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
         
            lij  = inddm(ind_loop(1) , j, k)

#ifdef _OPENMP4
CCCC!$OMP simd
#else
CDIR$ IVDEP
#endif
            do  l = lij, lij + ind_loop(2) - ind_loop(1)

              ropr2 = rop(l,2)*rot(1,1)+rop(l,3)*rot(1,2)+
     &                rop(l,4)*rot(1,3)
              ropr3 = rop(l,2)*rot(2,1)+rop(l,3)*rot(2,2)+
     &                rop(l,4)*rot(2,3)
              ropr4 = rop(l,2)*rot(3,1)+rop(l,3)*rot(3,2)+
     &                rop(l,4)*rot(3,3)
            
              rop(l,1) = rop(l,1)*c1
              rop(l,2) = ropr2*c1
              rop(l,3) = ropr3*c1
              rop(l,4) = ropr4*c1
              rop(l,5) = rop(l,5)*c1
            enddo

         enddo
         enddo

      else

         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)

            lij  = inddm(ind_loop(1) , j, k)

#ifdef _OPENMP4
CCCC!$OMP simd
#else
CDIR$ IVDEP
#endif
            do  l = lij, lij + ind_loop(2) - ind_loop(1)

              ropr2 = rop(l,2)*rot(1,1)+rop(l,3)*rot(1,2)+
     &                rop(l,4)*rot(1,3)
              ropr3 = rop(l,2)*rot(2,1)+rop(l,3)*rot(2,2)+
     &                rop(l,4)*rot(2,3)
              ropr4 = rop(l,2)*rot(3,1)+rop(l,3)*rot(3,2)+
     &                rop(l,4)*rot(3,3)
            
              rop(l,1) = rop(l,1)*c1
              rop(l,2) = ropr2*c1
              rop(l,3) = ropr3*c1
              rop(l,4) = ropr4*c1
              rop(l,5) = rop(l,5)*c1
            enddo

         enddo
         enddo

      endif

c      write(*,*)'idir=', idir,inc1,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop=', ind_loop

      end
