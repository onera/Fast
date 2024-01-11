c***********************************************************************
c     $Date: 2010-04-06 11:08:16 +0200 (Tue, 06 Apr 2010) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine nature_geom_dom(nijk_xyz,
     &                           ndimt_xyz, x,y,z, degener,
     &                           lale, param_real,
     &                           nijk_mtr, ndimdx_mtr,
     &                           neq_ij, neq_k, type_zone)
c***********************************************************************
c     AC 
c_A    determine si le domaine est 3D, 3D + k direction homoge, ou 2D
c
c     INP
c_I    ndom     : numero du domaine courant
c_I    ndimt_xyz    : dimension x,y,z
c_I    x,y,z    : coordonnees des centres de cellules
c
c     OUT
c_O    neq_ij   : nombre de composante pour les metriques des facettes i et j
c_O    neq_k    : nombre de composante pour les metriques des facettes k
c***********************************************************************
      implicit none

#include "FastC/param_solver.h"
      REAL_E souszero
c      parameter(souszero=1e-10)
      parameter(souszero=2.e-8)


      INTEGER_E ndimt_xyz,neq_ij,neq_k,type_zone,ndimdx_mtr,lale,
     &  nijk_xyz(5),nijk_mtr(5),degener(ndimt_xyz)

      REAL_E x(ndimt_xyz),y(ndimt_xyz),z(ndimt_xyz), param_real(0:*)

c Var loc
      logical xflag,yflag,xcart,xcte,ycart,ycte,zcart,zcte
      INTEGER_E ind0,ind,ind1,ind2,i,j,k,kcible,jmax,kmax,l

      REAL_E xtest,ytest,ztest1,ztest2,ztest3,ztest4,xm,ym,zm,cible,
     & degen1,degen2,degen3,degen4

#include "FastC/formule_xyz.h"

         !do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5) 
         !do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4) 
         !do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4) 
         ! ind  = indcg(i  , j  , k  )
         ! write(*,'(3f18.10,3i7)')x(ind),y(ind),z(ind),i,j,k
         !enddo
         !enddo
         !enddo

         !recherche maille ecrasee par application CL
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5) - 1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4) - 1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4) - 1


          !direction i
          ind  = indcg(i  , j+1  , k  )
          ind1 = indcg(i+1, j+1  , k  )
          degen2 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          ind  = indcg(i  , j+1  , k+1  )
          ind1 = indcg(i+1, j+1  , k+1  )
          degen3 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          ind  = indcg(i  , j    , k+1  )
          ind1 = indcg(i+1, j    , k+1  )
          degen4 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          ind  = indcg(i  , j  , k  )
          ind1 = indcg(i+1, j  , k  )
          degen1 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          if(degen1.le.souszero.and.degen2.le.souszero.and.
     &       degen3.le.souszero.and.degen4.le.souszero) then

               degener(ind) = 1
          else 
               degener(ind) = 0
          endif


          !direction j
          ind  = indcg(i+1, j    , k  )
          ind1 = indcg(i+1, j+1  , k  )
          degen2 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

!          if(degen2.le.souszero) write(*,'(a,6f18.10)')'degen2',
!     &              x(ind),x(ind1), y(ind),y(ind1),z(ind),z(ind1)

          ind  = indcg(i+1, j    , k+1  )
          ind1 = indcg(i+1, j+1  , k+1  )
          degen3 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

!          if(degen3.le.souszero) write(*,'(a,6f18.10)')'degen3',
!     &              x(ind),x(ind1), y(ind),y(ind1),z(ind),z(ind1)

          ind  = indcg(i  , j    , k+1  )
          ind1 = indcg(i  , j+1  , k+1  )
          degen4 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

!          if(degen4.le.souszero) write(*,'(a,6f18.10)')'degen4',
!     &              x(ind),x(ind1), y(ind),y(ind1),z(ind),z(ind1)
c     &              x(ind)-x(ind1), y(ind)-y(ind1),z(ind)-z(ind1)

          ind  = indcg(i  , j  , k  )
          ind1 = indcg(i  , j+1, k  )
          degen1 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

!          if(degen1.le.souszero) write(*,'(a,6f18.10)')'degen1',
!     &              x(ind),x(ind1), y(ind),y(ind1),z(ind),z(ind1)
c     &              x(ind)-x(ind1), y(ind)-y(ind1),z(ind)-z(ind1)

          if(degen1.le.souszero.and.degen2.le.souszero.and.
     &       degen3.le.souszero.and.degen4.le.souszero.and.
     &       degener(ind).eq.0) then

               degener(ind) = 1
               !write(*,*)'degen J', i,j,k
               !write(*,'(4f18.10)')degen1,degen2,degen3,degen4
          endif


          !direction k
          ind  = indcg(i+1, j    , k  )
          ind1 = indcg(i+1, j    , k+1)
          degen2 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          ind  = indcg(i+1, j+1  , k    )
          ind1 = indcg(i+1, j+1  , k+1  )
          degen3 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          ind  = indcg(i  , j+1  , k    )
          ind1 = indcg(i  , j+1  , k+1  )
          degen4 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          ind  = indcg(i  , j  , k  )
          ind1 = indcg(i  , j  , k+1)
          degen1 = sqrt(  (x(ind)-x(ind1))**2
     &                  + (y(ind)-y(ind1))**2 
     &                  + (z(ind)-z(ind1))**2  )

          if(degen1.le.souszero.and.degen2.le.souszero.and.
     &       degen3.le.souszero.and.degen4.le.souszero.and.
     &       degener(ind).eq.0) then

               degener(ind) = 1
               !write(*,*)'degen K', i,j,k
               !write(*,*)degen1,degen2,degen3,degen4
          endif

         enddo
         enddo
         enddo




      !si domaine 2d
      IF(nijk_xyz(3).eq.2) THEN

        type_zone = 3
        neq_ij    = 2
        neq_k     = 0

        nijk_mtr(1) =  1
        nijk_mtr(2) =  nijk_xyz(1)
        nijk_mtr(3) =  0
        nijk_mtr(4) =  nijk_xyz(4)
        nijk_mtr(5) =  nijk_xyz(5)

        ndimdx_mtr  = nijk_xyz(1)*nijk_xyz(2)  ! taille tableau metric = ni*nj

c       !Protection pour eviter deplacement axe z sur dom 2D
c       if(lale.and.
c     &   (omega(1).ne.0.or.omega(2).ne.0.or.vtrans(3).ne.0)) then
c
c        write(*,*)'Erreur: rotation axe x/y ou translation axe z 
c     & impossible sur dom 2D' 
c        call error('nature_geom_dom$',70,1)
c       endif

      ELSE !domaine 3D

       xcart    = .true.
       ycart    = .true.
       zcart    = .true.
       xcte     = .true.
       ycte     = .true.
       zcte     = .true.

       xflag = .true.
       yflag = .true.

       xm = 0.
       ym = 0.
       zm = 0.
 
       !!!
       !!!
       !!! direction k
       !!!
       !!!

         ! cartesienne ?
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5) - 1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4) - 1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4) - 1

          ind0 = indcg(1  , 1    , k  )
          ind1 = indcg(i  , j    , k  )
          ind2 = indcg(i+1, j+1  , k  )

          ztest1= abs( z(ind0)-z(ind1) )
          ztest2= abs( z(ind0)-z(ind2) )

          ind0 = indcg(1  , 1    , k+1  )
          ind1 = indcg(i  , j    , k+1  )
          ind2 = indcg(i+1, j+1  , k+1  )

          ztest3= abs( z(ind0)-z(ind1) )
          ztest4= abs( z(ind0)-z(ind2) )

          ind  = indcg(i  , j  , k  )

          if((ztest1.ge.souszero.or.ztest2.ge.souszero.or.
     &        ztest3.ge.souszero.or.ztest4.ge.souszero).and.
     &        degener(ind).eq.0) then
            zcart= .false.
          endif
         enddo
         enddo
         enddo
         ! pas constant  ?
         cible = z(indcg(1,1,2) )-z( indcg(1,1,1) )
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5)-1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4)-1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4)-1

          ind  = indcg(i,j,k)
          ind1 = indcg(i,j,k+1)

          ztest1= abs((z(ind1)-z(ind)) - cible)
          if(ztest1.ge.souszero.and.degener(ind).eq.0) then
            zcte = .false.
          endif
         enddo
         enddo
         enddo   

       !!!
       !!!
       !!! direction j
       !!!
       !!!

         ! cartesienne ?
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5)-1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4)-1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4)-1

          ind0 = indcg(1  ,j , 1  )
          ind1 = indcg(i  ,j , k  )
          ind2 = indcg(i+1,j , k+1)

          ztest1= abs( y(ind0)-y(ind1) )
          ztest2= abs( y(ind0)-y(ind2) )

          ind0 = indcg(1  ,j+1 , 1  )
          ind1 = indcg(i  ,j+1 , k  )
          ind2 = indcg(i+1,j+1 , k+1)

          ztest3= abs( y(ind0)-y(ind1) )
          ztest4= abs( y(ind0)-y(ind2) )

          ind  = indcg(i,j,k)

          if((ztest1.ge.souszero.or.ztest2.ge.souszero.or.
     &        ztest3.ge.souszero.or.ztest4.ge.souszero).and.
     &       degener(ind).eq.0) then
            ycart= .false.
          endif

         enddo
         enddo
         enddo
         ! pas constant  ?
         cible = y(indcg(1,2,1) )-y( indcg(1,1,1) )
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5)-1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4)-1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4)-1

          ind  = indcg(i,j  ,k)
          ind1 = indcg(i,j+1,k)

          ztest1= abs((y(ind1)-y(ind)) - cible)
          if(ztest1.ge.souszero.and.degener(ind).eq.0) then
            ycte = .false.
          endif
         enddo
         enddo
         enddo

       !!!
       !!!
       !!! direction i
       !!!
       !!!

         ! cartesienne ?
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5)-1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4)-1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4)-1

          ind0 = indcg(i  ,1   , 1  )
          ind1 = indcg(i  ,j   , k  )
          ind2 = indcg(i  ,j+1 , k+1)

          ztest1= abs( x(ind0)-x(ind1) )
          ztest2= abs( x(ind0)-x(ind2) )

          ind0 = indcg(i+1  ,1   , 1  )
          ind1 = indcg(i+1  ,j   , k  )
          ind2 = indcg(i+1  ,j+1 , k+1)

          ztest3= abs( x(ind0)-x(ind1) )
          ztest4= abs( x(ind0)-x(ind2) )

          ind  = indcg(i,j,k)
          if((ztest1.ge.souszero.or.ztest2.ge.souszero.or.
     &        ztest3.ge.souszero.or.ztest4.ge.souszero).and.
     &       degener(ind).eq.0) then
              xcart= .false.
          endif
         enddo
         enddo
         enddo
         ! pas constant  ?
         cible = x(indcg(2,1,1) )-x( indcg(1,1,1) )
         do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5)-1
         do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4)-1
         do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4)-1

          ind  = indcg(i  ,j  ,k)
          ind1 = indcg(i+1,j  ,k)

          ztest1= abs((x(ind1)-x(ind)) - cible)
          if(ztest1.ge.souszero.and.degener(ind).eq.0) then
            xcte = .false.
          endif
         enddo
         enddo
         enddo

          !verifie si mouvement de rotation 2d/3d ou translation
          if (lale.eq.1) then
            ztest2=  param_real(ROT_OMEGA  )* param_real(ROT_OMEGA  )
     &           + param_real(ROT_OMEGA+1)* param_real(ROT_OMEGA+1)

            ztest1= ztest2+param_real(ROT_OMEGA+2)*param_real(ROT_OMEGA+2)

            ztest2= sqrt(ztest2)
            ztest1= sqrt(ztest1)
          else
            ztest1 = 0.
            ztest2 = 0.
          endif

         ! domaine cartesien pas constant
         if(xcart.and.xcte.and.ycart.and.ycte
     &                    .and.zcart.and.zcte
     &      .and.lale.le.1.and.ztest1.le.1.e-7) then

                   type_zone = 2
                   neq_ij    = 1
                   neq_k     = 1

                   nijk_mtr(1) =  0
                   nijk_mtr(2) =  0
                   nijk_mtr(3) =  0
                   nijk_mtr(4) =  0
                   nijk_mtr(5) =  0
                   ndimdx_mtr  =  1
  
         elseif(zcart.and.zcte) then

             do k=1-nijk_xyz(5), nijk_xyz(3) -nijk_xyz(5)-1
             do j=1-nijk_xyz(4), nijk_xyz(2) -nijk_xyz(4)-1
             do i=1-nijk_xyz(4), nijk_xyz(1) -nijk_xyz(4)-1

               ind0 = indcg(i  ,j  ,1)
               ind1 = indcg(i  ,j  ,k)
               ztest1= abs( x(ind0)-x(ind1) )

               ind0 = indcg(i+1  ,j+1  ,1)
               ind1 = indcg(i+1  ,j+1  ,k)
               ztest2= abs( x(ind0)-x(ind1) )

               ind0 = indcg(i  ,j  , 1  )
               ind1 = indcg(i  ,j  , k+1)
               ztest3= abs( x(ind0)-x(ind1) )

               ind0 = indcg(i+1  ,j+1  , 1  )
               ind1 = indcg(i+1  ,j+1  , k+1)
               ztest4= abs( x(ind0)-x(ind1) )
               
               ind  = indcg(i  ,j  ,k)

               if((ztest1.ge.souszero.or.ztest2.ge.souszero.or.
     &             ztest3.ge.souszero.or.ztest4.ge.souszero).and.
     &              degener(ind).eq.0) then
                   xflag= .false.
               endif

               ind0 = indcg(i  ,j  ,1)
               ind1 = indcg(i  ,j  ,k)
               ztest1= abs( y(ind0)-y(ind1) )

               ind0 = indcg(i+1  ,j+1  ,1)
               ind1 = indcg(i+1  ,j+1  ,k)
               ztest2= abs( y(ind0)-y(ind1) )

               ind0 = indcg(i  ,j  , 1  )
               ind1 = indcg(i  ,j  , k+1)
               ztest3= abs( y(ind0)-y(ind1) )

               ind0 = indcg(i+1  ,j+1  , 1  )
               ind1 = indcg(i+1  ,j+1  , k+1)
               ztest4= abs( y(ind0)-y(ind1) )

               if((ztest1.ge.souszero.or.ztest2.ge.souszero.or.
     &             ztest3.ge.souszero.or.ztest4.ge.souszero).and.
     &              degener(ind).eq.0) then
                   yflag= .false.
               endif

             enddo
             enddo
             enddo

             if(xflag.and.yflag.and.lale.le.1.and.ztest2.le.1.e-7) then
                 !domaine  3d homogene avec rotation eventuelle sur axe z
                 type_zone = 1
                 neq_ij    = 2
                 neq_k     = 1

                 nijk_mtr(1) =  1
                 nijk_mtr(2) =  nijk_xyz(1)
                 nijk_mtr(3) =  0
                 nijk_mtr(4) =  nijk_xyz(4)
                 nijk_mtr(5) =  0
                 ndimdx_mtr  =  nijk_xyz(1)*nijk_xyz(2)
             else
                 !domaine 3d general
                 type_zone = 0
                 neq_ij    = 3
                 neq_k     = 3

                 nijk_mtr(1) =  1
                 nijk_mtr(2) =  nijk_xyz(1)
                 nijk_mtr(3) =  nijk_xyz(1)*nijk_xyz(2)
                 nijk_mtr(4) =  nijk_xyz(4)
                 nijk_mtr(5) =  nijk_xyz(5)
                 ndimdx_mtr  =  nijk_xyz(1)*nijk_xyz(2)*nijk_xyz(3)
             endif

         else    !domaine  3d general
                 type_zone = 0
                 neq_ij    = 3
                 neq_k     = 3

                 nijk_mtr(1) =  1
                 nijk_mtr(2) =  nijk_xyz(1)
                 nijk_mtr(3) =  nijk_xyz(1)*nijk_xyz(2)
                 nijk_mtr(4) =  nijk_xyz(4)
                 nijk_mtr(5) =  nijk_xyz(5)
                 ndimdx_mtr  =  nijk_xyz(1)*nijk_xyz(2)*nijk_xyz(3)
         endif
         
c       IF(lale.and.(omega(1).ne.0.or.omega(2).ne.0                   
c    &                           .or.vtrans(3).ne.0)) xflag=.false. !pas d'optimisation si ale axe z

       !xflag=.false.

       !write(*,*)'flag   ',xflag,yflag,zflag,ndom
       !write(*,*)'flag er',xm,ym,zm

      ENDIF!domaine 3D/2D
     
      end
