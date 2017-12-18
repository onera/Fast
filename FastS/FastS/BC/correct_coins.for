c***********************************************************************
c     $Date: 2012-02-21 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: AlferezNicolas $
c***********************************************************************
      subroutine correct_coins(ndom,param_int, ind_loop, rop)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"


      INTEGER_E ndom, param_int(0:*),ind_loop(6)
      REAL_E rop( param_int(NDIMDX), param_int(NEQ))

c Var Local
      INTEGER_E  i,j,k,l1,l11,l12,l2,l13,ne, ind_loop_loc(6),
     & l,i1,j1,k1,i0,j0,k0,inddm_loc,itest(4),is,idir
      REAL_E c3

#include "FastS/formule_param.h"

        ind_loop_loc(1)=max(1,ind_loop(1))
        ind_loop_loc(3)=max(1,ind_loop(3))
        ind_loop_loc(5)=max(1,ind_loop(5))
        ind_loop_loc(2)=min(param_int(IJKV),ind_loop(2))
        ind_loop_loc(4)=min(param_int(IJKV+1),ind_loop(4))
        ind_loop_loc(6)=min(param_int(IJKV+2),ind_loop(6))
        !face I = Imin et Imax: arrete J et K
        itest(1)=1
        itest(2)=param_int(IJKV)

        do idir=1,2 

         is = 1
         if(idir.eq.1) is = -1

         if(ind_loop_loc(idir).eq.itest(idir)) then

           if(ind_loop_loc(3).eq.1) then !Bande jmin

            do ne = 1, param_int(NEQ)
            do k  = ind_loop_loc(5),ind_loop_loc(6)
            do j  = 1,param_int(NIJK+3)
            do i  = 1,param_int(NIJK+3)
 
             i1 = min(param_int(NIJK+3),i+1)
             j1 = min(param_int(NIJK+3),j+1)

             l1 =inddm( ind_loop_loc(idir)+i*is , ind_loop_loc(3)-j, k)
             l11=inddm( ind_loop_loc(idir)+i1*is, ind_loop_loc(3)  , k)
             l12=inddm( ind_loop_loc(idir)      , ind_loop_loc(3)-j1,k)

             rop(l1,ne)=(rop(l11,ne)+rop(l12,ne))*0.5

            enddo
            enddo
            enddo
            enddo
           endif
           if(ind_loop_loc(4).eq.param_int(IJKV+1)) then !Bande jmax

            do ne = 1, param_int(NEQ)
            do k=ind_loop_loc(5),ind_loop_loc(6)
            do j=1,param_int(NIJK+3)
            do i=1,param_int(NIJK+3)
 
             i1 = min(param_int(NIJK+3),i+1)
             j1 = min(param_int(NIJK+3),j+1)

             l1 =inddm( ind_loop_loc(idir)+i*is , ind_loop_loc(4)+j, k)
             l11=inddm( ind_loop_loc(idir)+i1*is, ind_loop_loc(4)  , k)
             l12=inddm( ind_loop_loc(idir)      , ind_loop_loc(4)+j1,k)

             rop(l1,ne)=(rop(l11,ne)+rop(l12,ne))*0.5
            enddo
            enddo
            enddo
            enddo
           endif

           if(ind_loop_loc(5).eq.1) then !Bande kmin

            do ne = 1, param_int(NEQ)
            do k=1,param_int(NIJK+4)
            do j=ind_loop_loc(3),ind_loop_loc(4)
            do i=1,param_int(NIJK+3)

             k1 = min(param_int(NIJK+4) , k+1)
             i1 = min(param_int(NIJK+3), i+1)

             l1 =inddm( ind_loop_loc(idir)+i*is ,j, ind_loop_loc(5)-k )
             l11=inddm( ind_loop_loc(idir)      ,j, ind_loop_loc(5)-k1)
             l12=inddm( ind_loop_loc(idir)+i1*is,j, ind_loop_loc(5)   )

             rop(l1,ne)=(rop(l11,ne)+rop(l12,ne))*0.5
            enddo
            enddo
            enddo
            enddo
           endif
           if(ind_loop_loc(6).eq.param_int(IJKV+2)) then !Bande kmax

            do ne = 1, param_int(NEQ)
            do k=1,param_int(NIJK+4)
            do j=ind_loop_loc(3),ind_loop_loc(4)
            do i=1,param_int(NIJK+3)

             k1 = min(param_int(NIJK+4) , k+1)
             i1 = min(param_int(NIJK+3), i+1)
                                   
             l1 =inddm( ind_loop_loc(idir)+i*is ,j, ind_loop_loc(6)+k )
             l11=inddm( ind_loop_loc(idir)      ,j, ind_loop_loc(6)+k1)
             l12=inddm( ind_loop_loc(idir)+i1*is,j, ind_loop_loc(6)   )

             rop(l1,ne)=(rop(l11,ne)+rop(l12,ne))*0.5
            enddo
            enddo
            enddo
            enddo
           endif
          endif! le sous domaine touche la Face Imin ou Imax du domaine
        enddo !idir Face i



         !face I = Imin et Imax: arrete J et K
        itest(3)=1
        itest(4)=param_int(IJKV+1)

        do idir=3,4 

         is = 1
         if(idir.eq.3) is = -1

         if(ind_loop_loc(idir).eq.itest(idir)) then

           if(ind_loop_loc(5).eq.1) then !Bande kmin

            do ne = 1, param_int(NEQ)
            do k=1,param_int(NIJK+4)
            do j=1,param_int(NIJK+3)
CDIR$ IVDEP
            do i=ind_loop_loc(1),ind_loop_loc(2)

             k1 = min(param_int(NIJK+4) , k+1)
             j1 = min(param_int(NIJK+3), j+1)

             l1 =inddm(i, ind_loop_loc(idir)+j*is , ind_loop_loc(5)-k )
             l11=inddm(i, ind_loop_loc(idir)      , ind_loop_loc(5)-k1)
             l12=inddm(i, ind_loop_loc(idir)+j1*is, ind_loop_loc(5)   )

             rop(l1,ne)=(rop(l11,ne)+rop(l12,ne))*0.5
            enddo
            enddo
            enddo
            enddo
           endif
           if(ind_loop_loc(6).eq.param_int(IJKV+2)) then !Bande kmax

            do ne = 1, param_int(NEQ)
            do k=1,param_int(NIJK+4)
            do j=1,param_int(NIJK+3)
CDIR$ IVDEP
            do i=ind_loop_loc(1),ind_loop_loc(2)

             k1 = min(param_int(NIJK+4) , k+1)
             j1 = min(param_int(NIJK+3), j+1)

             l1 =inddm(i, ind_loop_loc(idir)+j*is , ind_loop_loc(6)+k )
             l11=inddm(i, ind_loop_loc(idir)      , ind_loop_loc(6)+k1)
             l12=inddm(i, ind_loop_loc(idir)+j1*is, ind_loop_loc(6)   )

             rop(l1,ne)=(rop(l11,ne)+rop(l12,ne))*0.5
            enddo
            enddo
            enddo
            enddo
           endif          
          endif! le sous domaine touche la Face Jmin ou Jmax du domaine
        enddo !idir Face J                         
 

         c3=1./3.

         !coin(1,1,1)
         i0 = ind_loop_loc(1)
         j0 = ind_loop_loc(3)
         k0 = ind_loop_loc(5)
         if(i0.eq.1.and.j0.eq.1.and.k0.eq.1) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)

            l  = inddm( i0-i , j0-j  , k0-k  )
            l11= inddm( i0   , j0    , k0-k1 )
            l12= inddm( i0-i1, j0    , k0    )
            l13= inddm( i0   , j0-j1 , k0    )
                   
            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

         !coin(iv,1,1)
         i0 = ind_loop_loc(2)
         j0 = ind_loop_loc(3)
         k0 = ind_loop_loc(5)
         if(i0.eq.param_int(IJKV).and.j0.eq.1.and.k0.eq.1) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)

            l  = inddm( i0+i , j0-j  , k0-k )
            l11= inddm( i0   , j0    , k0-k1)
            l12= inddm( i0+i1, j0    , k0   )
            l13= inddm( i0   , j0-j1 , k0   )

            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

        !coin(1,jv,1)
         i0 = ind_loop_loc(1)
         j0 = ind_loop_loc(4)
         k0 = ind_loop_loc(5)
         if(i0.eq.1.and.j0.eq.param_int(IJKV+1).and.k0.eq.1) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)

            l  = inddm( i0-i , j0+j  , k0-k )
            l11= inddm( i0   , j0    , k0-k1)
            l12= inddm( i0-i1, j0    , k0   )
            l13= inddm( i0   , j0+j1 , k0   )

            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

         !coin(iv,jv,1)
         i0 = ind_loop_loc(2)
         j0 = ind_loop_loc(4)
         k0 = ind_loop_loc(5)
         if(i0.eq.param_int(IJKV).and.j0.eq.param_int(IJKV+1)
     &                           .and.k0.eq.1) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)
            
            l  = inddm( i0+i , j0+j  , k0-k )
            l11= inddm( i0   , j0    , k0-k1)
            l12= inddm( i0+i1, j0    , k0   )
            l13= inddm( i0   , j0+j1 , k0   )
                   
            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

         !coin(1,1,kv)
         i0 = ind_loop_loc(1)
         j0 = ind_loop_loc(3)
         k0 = ind_loop_loc(6)
         if(i0.eq.1.and.j0.eq.1.and.k0.eq.param_int(IJKV+2)) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)

            l  = inddm( i0-i , j0-j  , k0+k  )
            l11= inddm( i0   , j0    , k0+k1 )
            l12= inddm( i0-i1, j0    , k0    )
            l13= inddm( i0   , j0-j1 , k0    )
                   
            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

         !coin(iv,1,kv)
         i0 = ind_loop_loc(2)
         j0 = ind_loop_loc(3)
         k0 = ind_loop_loc(6)
         if(i0.eq.param_int(IJKV).and.j0.eq.1
     &                           .and.k0.eq.param_int(IJKV+2)) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)

            l  = inddm( i0+i , j0-j  , k0+k )
            l11= inddm( i0   , j0    , k0+k1)
            l12= inddm( i0+i1, j0    , k0   )
            l13= inddm( i0   , j0-j1 , k0   )

            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

        !coin(1,jv,kv)
         i0 = ind_loop_loc(1)
         j0 = ind_loop_loc(4)
         k0 = ind_loop_loc(6)
         if(i0.eq.1.and.j0.eq.param_int(IJKV+1)
     &             .and.k0.eq.param_int(IJKV+2)) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)

            l  = inddm( i0-i , j0+j  , k0+k )
            l11= inddm( i0   , j0    , k0+k1)
            l12= inddm( i0-i1, j0    , k0   )
            l13= inddm( i0   , j0+j1 , k0   )

            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

         !coin(iv,jv,kv)
         i0 = ind_loop_loc(2)
         j0 = ind_loop_loc(4)
         k0 = ind_loop_loc(6)
         if(i0.eq.param_int(IJKV).and.j0.eq.param_int(IJKV+1)
     &                           .and.k0.eq.param_int(IJKV+2)) then
           do ne = 1, param_int(NEQ)
           do k=1,param_int(NIJK+4)
           do j=1,param_int(NIJK+3)
           do i=1,param_int(NIJK+3)

            i1 = min(param_int(NIJK+3), i+1)
            j1 = min(param_int(NIJK+3), j+1)
            k1 = min(param_int(NIJK+4) , k+1)
            
            l  = inddm( i0+i , j0+j  , k0+k )
            l11= inddm( i0   , j0    , k0+k1)
            l12= inddm( i0+i1, j0    , k0   )
            l13= inddm( i0   , j0+j1 , k0   )
                   
            rop(l,ne)=(rop(l11,ne)+rop(l12,ne)+rop(l13,ne))*c3
           enddo
           enddo
           enddo
           enddo
         endif

      end
