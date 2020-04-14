c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine dist_extrap(ndimdx,ndimt_xyz, nijk, nijk_xyz, ind_dm,
     &                      degener, dist)
c***********************************************************************
      implicit none

      INTEGER_E ndimdx,ndimt_xyz,
     &   ind_dm(6),nijk(5),nijk_xyz(5), degener(ndimt_xyz)

      REAL_E dist(ndimdx)

C_LOCAL
      INTEGER_E npass,i,j,k,imin,jmin,kmin,imax,jmax,kmax,lmtr,lmtr0,
     & lmtri,lmtrj,lmtr0i,lmtr0j,lmtrk,lmtr0k,ind_loop(6),ne,ind0,ind1,
     & ind3

#include "FastC/formule.h"
#include "FastC/formule_xyz.h"


      do npass=1,2

      if(npass.eq.1) then
        ind_loop(:)=ind_dm(:)
      else
       ind_loop(1)= ind_dm(1)-nijk(4)
       ind_loop(2)= ind_dm(2)+nijk(4)
       ind_loop(3)= ind_dm(3)-nijk(4)
       ind_loop(4)= ind_dm(4)+nijk(4)
       ind_loop(5)= ind_dm(5)-nijk(5)
       ind_loop(6)= ind_dm(6)+nijk(5)
      endif

      if(nijk(5).ne.0) then

         !!
         !!
         !!Kmin
         !!
         !!
         do k= 0        , ind_dm(5)-nijk(5), -1
         do j= ind_loop(3), ind_loop(4)
         do i= ind_loop(1), ind_loop(2)

             ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr =  inddm( i, j, k)
               lmtr0=  inddm( i, j, 1-k)

               dist(lmtr) = dist(lmtr0)
            endif
         enddo
         enddo
         enddo

         !!
         !!
         !!Kmax
         !!
         !!
         do k= ind_dm(6)+1          , ind_dm(6)+nijk(5) 
         do j= ind_loop(3), ind_loop(4)
         do i= ind_loop(1), ind_loop(2)

             ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr =  inddm( i, j, k                 )
               lmtr0=  inddm( i, j, 2*ind_dm(6) +1 -k )

               dist(lmtr) = dist(lmtr0)
            endif
         enddo
         enddo
         enddo


      endif

       !!
       !!
       !  Correction metrique plan JMIN fictif
       !!
       !!
         do k= ind_loop(5), ind_loop(6)
         do j= 0        , ind_dm(3)-nijk(4), -1
         do i= ind_loop(1), ind_loop(2)

            ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr =  inddm( i, j  ,k)
               lmtr0=  inddm( i, 1-j,k)

               dist(lmtr) = dist(lmtr0)
            endif
         enddo
         enddo
         enddo

       !!
       !!
       !  Correction metrique plan JMAX fictif
       !!
       !!
         do k= ind_loop(5), ind_loop(6)
         do j= ind_dm(4)+1          , ind_dm(4)+nijk(4)
         do i= ind_loop(1), ind_loop(2)

             ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr  =  inddm( i  , j                 , k)
               lmtr0 =  inddm( i  , 2*ind_dm(4) +1 -j , k)

               dist(lmtr) = dist(lmtr0)
            endif
         enddo
         enddo
         enddo

       !!
       !!
       !  Correction metrique plan IMIN fictif
       !!
       !!
         do k= ind_loop(5), ind_loop(6)
         do j= ind_loop(3), ind_loop(4)
         do i= 0        , ind_dm(1)-nijk(4), -1

             ind0 = indcg(i, j , k)

             if(   degener(ind0).eq.1 )  then

               lmtr =  inddm( i  ,j,k)
               lmtr0=  inddm( 1-i,j,k)

               dist(lmtr) = dist(lmtr0)
            endif
         enddo
         enddo
         enddo

       !!
       !!
       !  Correction metrique plan IMAX fictif
       !!
       !!
         do k= ind_loop(5), ind_loop(6)
         do j= ind_loop(3), ind_loop(4)
         do i= ind_dm(2)+1 , ind_dm(2)+nijk(4)

             ind0 = indcg(i, j , k)

             if(   degener(ind0).eq.1 )  then

               lmtr  =  inddm( i                 , j, k)
               lmtr0 =  inddm( 2*ind_dm(2)+ 1 -i , j, k)

               dist(lmtr) = dist(lmtr0)
            endif
         enddo
         enddo
         enddo

      enddo !npass

      end

