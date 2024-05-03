c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine tijk_extrap(ndimdx_mtr,ndimt_xyz,
     &                      nijk_xyz, nijk_mtr,
     &                      neq_ij,neq_k,
     &                      ind_dm,
     &                      degener,
     &                      ti,tj,tk, ti0,tj0,tk0, vol)
c***********************************************************************
      implicit none

      INTEGER_E ndimdx_mtr,ndimt_xyz,neq_ij,neq_k,
     &   ind_dm(6),nijk_xyz(5),nijk_mtr(5), degener(ndimt_xyz)

      REAL_E ti(ndimdx_mtr,neq_ij), ti0(ndimdx_mtr,neq_ij)   ! metric
      REAL_E tj(ndimdx_mtr,neq_ij), tj0(ndimdx_mtr,neq_ij)
      REAL_E tk(ndimdx_mtr,neq_k ), tk0(ndimdx_mtr,neq_k )
      REAL_E vol(ndimdx_mtr)

C_LOCAL
      INTEGER_E npass,i,j,k,imin,jmin,kmin,imax,jmax,kmax,lmtr,lmtr0,
     & lmtri,lmtrj,lmtr0i,lmtr0j,lmtrk,lmtr0k,ind_loop(6),ne,ind0,ind1,
     & ind3

#include "FastC/formule_mtr.h"
#include "FastC/formule_xyz.h"


      do npass=1,2

      if(npass.eq.1) then
        ind_loop(:)=ind_dm(:)
      else
       ind_loop(1)= ind_dm(1)-nijk_mtr(4)
       ind_loop(2)= ind_dm(2)+nijk_mtr(4)
       ind_loop(3)= ind_dm(3)-nijk_mtr(4)
       ind_loop(4)= ind_dm(4)+nijk_mtr(4)
       ind_loop(5)= ind_dm(5)-nijk_mtr(5)
       ind_loop(6)= ind_dm(6)+nijk_mtr(5)
      endif

      if(nijk_mtr(5).ne.0) then

         !!
         !!
         !!Kmin
         !!
         !!
         do k= 0        , ind_dm(5)-nijk_mtr(5), -1
         do j= ind_loop(3), ind_loop(4)
         do i= ind_loop(1), ind_loop(2)

             imin = max( ind_dm(1)-nijk_mtr(4) , i-1 )
             jmin = max( ind_dm(3)-nijk_mtr(4) , j-1 )
             imax = min (ind_dm(2)+nijk_mtr(4) , i+1 )
             jmax = min (ind_dm(4)+nijk_mtr(4) , j+1 )

             ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr =  indmtr( i, j, k)
               lmtr0=  indmtr( i, j, 1-k)

               ind1 = indcg(imin, j ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   ti( lmtr , ne) = ti( lmtr0 , ne)
                   ti0(lmtr , ne) = ti0(lmtr0 , ne)
                 enddo
               endif
               ind1 = indcg(i, jmin ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   tj( lmtr , ne) = tj( lmtr0 , ne)
                   tj0(lmtr , ne) = tj0(lmtr0 , ne)
                 enddo
               endif

               lmtri =  indmtr( i+1, j, k)
               lmtr0i=  indmtr( i+1, j, 1-k)
               ind1  =  indcg(imax,j,k)

               if(degener(ind1).eq.1) then

                 do ne =1, neq_ij
                   ti( lmtri, ne) = ti( lmtr0i, ne)
                   ti0(lmtri, ne) = ti0(lmtr0i, ne)
                 enddo
               endif

               lmtrj =  indmtr( i  , j+1, k   )
               lmtr0j=  indmtr( i  , j+1, 1-k )
               ind1  =  indcg(i, jmax ,k  )

               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   tj( lmtrj , ne) = tj( lmtr0j , ne)
                   tj0(lmtrj , ne) = tj0(lmtr0j , ne)
                 enddo
               endif

               !!Reconstruction tk pour assurer conservation geom
               lmtrk =  indmtr( i  ,j, k+1)
               if(neq_k.ne.3) then

                 do ne =1, neq_k
                 tk( lmtr, ne) = tk( lmtrk,ne)
                 tk0(lmtr, ne) = tk0(lmtrk,ne)
                 enddo

               else
                 do ne =1, neq_ij
                 tk( lmtr, ne)= tk(lmtrk,ne)- ti(lmtr,ne)+ ti(lmtri,ne)
     &                                      - tj(lmtr,ne)+ tj(lmtrj,ne)
                 tk0(lmtr, ne)=tk0(lmtrk,ne)-ti0(lmtr,ne)+ ti0(lmtri,ne)
     &                                      -tj0(lmtr,ne)+ tj0(lmtrj,ne)
                 enddo

              endif

               vol(lmtr) = vol(lmtr0)

             endif
         enddo
         enddo
         enddo

         !!
         !!
         !!Kmax
         !!
         !!
         do k= ind_dm(6)+1          , ind_dm(6)+nijk_mtr(5) 
         do j= ind_loop(3), ind_loop(4)
         do i= ind_loop(1), ind_loop(2)

             imin = max( ind_dm(1)-nijk_mtr(4) , i-1 )
             jmin = max( ind_dm(3)-nijk_mtr(4) , j-1 )
             imax = min (ind_dm(2)+nijk_mtr(4) , i+1 )
             jmax = min (ind_dm(4)+nijk_mtr(4) , j+1 )

             ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr =  indmtr( i, j, k                 )
               lmtr0=  indmtr( i, j, 2*ind_dm(6) +1 -k )

               ind1 = indcg(imin, j ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   ti( lmtr , ne) = ti( lmtr0 , ne)
                   ti0(lmtr , ne) = ti0(lmtr0 , ne)
                 enddo
               endif
               ind1 = indcg(i, jmin ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   tj( lmtr , ne) = tj( lmtr0 , ne)
                   tj0(lmtr , ne) = tj0(lmtr0 , ne)
                 enddo
               endif

               lmtri =  indmtr( i+1, j, k                )
               lmtr0i=  indmtr( i+1, j, 2*ind_dm(6) +1 -k)
               ind1  =  indcg(imax,j,k)

               if(degener(ind1).eq.1) then

                 do ne =1, neq_ij
                   ti( lmtri, ne) = ti( lmtr0i, ne)
                   ti0(lmtri, ne) = ti0(lmtr0i, ne)
                 enddo
               endif

               lmtrj =  indmtr( i  , j+1, k                 )
               lmtr0j=  indmtr( i  , j+1, 2*ind_dm(6) +1 -k )
               ind1  =  indcg(i, jmax ,k  )

               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   tj( lmtrj , ne) = tj( lmtr0j , ne)
                   tj0(lmtrj , ne) = tj0(lmtr0j , ne)
                 enddo
               endif

               !!Reconstruction ti pour assurer conservation geom
               lmtrk =  indmtr( i  ,j, k+1)
               if(neq_k.ne.3) then

                 do ne =1, neq_k
                 tk( lmtrk, ne)=tk( lmtr,ne)
                 tk0(lmtrk, ne)=tk0(lmtr,ne)
                 enddo

               else
                 do ne =1, neq_ij
                 tk( lmtrk, ne)= tk(lmtr,ne)+ ti(lmtr,ne)- ti(lmtri,ne)
     &                                      + tj(lmtr,ne)- tj(lmtrj,ne)
                 tk0(lmtrk, ne)=tk0(lmtr,ne)+ti0(lmtr,ne)- ti0(lmtri,ne)
     &                                      +tj0(lmtr,ne)- tj0(lmtrj,ne)
                 enddo
               endif

               vol(lmtr) = vol(lmtr0)

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
         do j= 0        , ind_dm(3)-nijk_mtr(4), -1
         do i= ind_loop(1), ind_loop(2)

             imin = max( ind_dm(1)-nijk_mtr(4) , i-1 )
             kmin = max( ind_dm(5)-nijk_mtr(5) , k-1 )
             imax = min (ind_dm(2)+nijk_mtr(4) , i+1 )
             kmax = min (ind_dm(6)+nijk_mtr(5) , k+1 )

             ind0 = indcg(i   ,j  ,k)

            if( degener(ind0).eq.1) then

               lmtr =  indmtr( i, j  ,k)
               lmtr0=  indmtr( i, 1-j,k)

               ind1 = indcg(imin, j ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   ti( lmtr , ne) = ti( lmtr0 , ne)
                   ti0(lmtr , ne) = ti0(lmtr0 , ne)
                 enddo
               endif

               ind1 = indcg(i, j    ,kmin)
               if(degener(ind1).eq.1) then
                 do ne =1, neq_k
                   tk (lmtr  , ne)  = tk ( lmtr0 , ne)
                   tk0(lmtr  , ne)  = tk0( lmtr0, ne)
                 enddo
               endif

               lmtri =  indmtr( i+1, j  , k)
               lmtr0i=  indmtr( i+1, 1-j, k)
               ind1  =  indcg(imax,j,k)

               if(degener(ind1).eq.1) then

                 do ne =1, neq_ij
                   ti( lmtri, ne) = ti( lmtr0i, ne)
                   ti0(lmtri, ne) = ti0(lmtr0i, ne)
                 enddo
               endif

               lmtrk =  indmtr( i  ,j   ,k+1)
               lmtr0k=  indmtr( i  ,1-j ,k+1)
               ind1  =  indcg(i,j,kmax)

               if(degener(ind1).eq.1) then

                 do ne =1, neq_k
                   tk(lmtrk , ne)  = tk( lmtr0k, ne)
                   tk0(lmtrk, ne)  = tk0(lmtr0k, ne)
                 enddo
               endif

               !!Reconstruction ti pour assurer conservation geom
               lmtrj =  indmtr( i  ,j+1,k)
               if(neq_k.ne.3) then

                 do ne =1, neq_ij
                 tj( lmtr, ne)=tj( lmtrj,ne)-  ti(lmtr,ne)+ ti(lmtri,ne)
                 tj0(lmtr, ne)=tj0(lmtrj,ne)- ti0(lmtr,ne)+ti0(lmtri,ne)
                 enddo

               else
                 do ne =1, neq_ij
                 tj( lmtr, ne)= tj(lmtrj,ne)- ti(lmtr,ne)+ ti(lmtri,ne)
     &                                      - tk(lmtr,ne)+ tk(lmtrk,ne)
                 tj0(lmtr, ne)=tj0(lmtrj,ne)-ti0(lmtr,ne)+ ti0(lmtri,ne)
     &                                      -tk0(lmtr,ne)+ tk0(lmtrk,ne)
                 enddo

              endif


               vol(lmtr) = vol(lmtr0)

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
         do j= ind_dm(4)+1          , ind_dm(4)+nijk_mtr(4)
         do i= ind_loop(1), ind_loop(2)

             imin = max( ind_dm(1)-nijk_mtr(4) , i-1 )
             kmin = max( ind_dm(5)-nijk_mtr(5) , k-1 )
             imax = min (ind_dm(2)+nijk_mtr(4) , i+1 )
             kmax = min (ind_dm(6)+nijk_mtr(5) , k+1 )

             ind0 = indcg(i   ,j  ,k)
             ind1 = indcg(i, j    ,kmin)
             ind3 = indcg(i,j,kmax)

            if( degener(ind0).eq.1) then

               lmtr  =  indmtr( i  , j                 , k)
               lmtr0 =  indmtr( i  , 2*ind_dm(4) +1 -j , k)

               ind1 = indcg(imin, j ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   ti( lmtr , ne) = ti( lmtr0 , ne)
                   ti0(lmtr , ne) = ti0(lmtr0 , ne)
                 enddo
               endif

               ind1 = indcg(i, j    ,kmin)
               if(degener(ind1).eq.1) then
                 do ne =1, neq_k
                   tk (lmtr  , ne)  = tk ( lmtr0 , ne)
                   tk0(lmtr  , ne)  = tk0( lmtr0, ne)
                 enddo
               endif

               lmtri =  indmtr( i+1, j                 , k)
               lmtr0i=  indmtr( i+1, 2*ind_dm(4) +1 -j , k)

               ind1 = indcg(imax,j,k)
               if(degener(ind1).eq.1) then

                 do ne =1, neq_ij
                   ti( lmtri, ne) = ti( lmtr0i, ne)
                   ti0(lmtri, ne) = ti0(lmtr0i, ne)
                 enddo
               endif

               lmtrk =  indmtr( i, j                 , k+1)
               lmtr0k=  indmtr( i, 2*ind_dm(4) +1 -j , k+1)

               ind1 = indcg(i,j,kmax)
               if(degener(ind1).eq.1) then

                 do ne =1, neq_k
                   tk(lmtrk , ne)  = tk( lmtr0k, ne)
                   tk0(lmtrk, ne)  = tk0(lmtr0k, ne)
                 enddo
               endif

               !!Reconstruction ti pour assurer conservation geom
               lmtrj  =  indmtr( i , j+1         , k)
               if(neq_k.ne.3) then

                 do ne =1, neq_ij
                 tj( lmtrj, ne)=tj( lmtr,ne)+  ti(lmtr,ne)- ti(lmtri,ne)
                 tj0(lmtrj, ne)=tj0(lmtr,ne)+ ti0(lmtr,ne)-ti0(lmtri,ne)
                 enddo

               else
                 do ne =1, neq_ij
                 tj( lmtrj, ne)= tj(lmtr,ne)+ ti(lmtr,ne)- ti(lmtri,ne)
     &                                      + tk(lmtr,ne)- tk(lmtrk,ne)
                 tj0(lmtrj, ne)=tj0(lmtr,ne)+ti0(lmtr,ne)- ti0(lmtri,ne)
     &                                      +tk0(lmtr,ne)- tk0(lmtrk,ne)
                 enddo
              endif

               vol(lmtr) = vol(lmtr0)

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
         do i= 0        , ind_dm(1)-nijk_mtr(4), -1


             jmin = max( ind_dm(3)-nijk_mtr(4) , j-1 )
             kmin = max( ind_dm(5)-nijk_mtr(5) , k-1 )
             jmax = min (ind_dm(4)+nijk_mtr(4) , j+1 )
             kmax = min (ind_dm(6)+nijk_mtr(5) , k+1 )

             ind0 = indcg(i, j , k)

             if(   degener(ind0).eq.1 )  then

               lmtr =  indmtr( i  ,j,k)
               lmtr0=  indmtr( 1-i,j,k)

               ind1 = indcg(i, jmin ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   tj( lmtr , ne) = tj( lmtr0 , ne)
                   tj0(lmtr , ne) = tj0(lmtr0 , ne)
                 enddo
               endif

               ind1 = indcg(i, j    ,kmin)
               if(degener(ind1).eq.1) then
                 do ne =1, neq_k
                   tk (lmtr  , ne)  = tk ( lmtr0 , ne)
                   tk0(lmtr  , ne)  = tk0( lmtr0, ne)
                 enddo
               endif

               lmtrj =  indmtr( i  ,j+1 ,k)
               lmtr0j=  indmtr( 1-i,j+1 ,k)

               ind1 = indcg(i,jmax,k)
               if(degener(ind1).eq.1) then

                 do ne =1, neq_ij
                   tj( lmtrj, ne) = tj( lmtr0j, ne)
                   tj0(lmtrj, ne) = tj0(lmtr0j, ne)
                 enddo
               endif

               lmtrk =  indmtr( i  ,j ,k+1)
               lmtr0k=  indmtr( 1-i,j ,k+1)

               ind1 = indcg(i,j,kmax)
               if(degener(ind1).eq.1) then

                 do ne =1, neq_k
                   tk( lmtrk, ne)  = tk( lmtr0k, ne)
                   tk0(lmtrk, ne)  = tk0(lmtr0k, ne)
                 enddo
               endif

               !!Reconstruction ti pour assurer conservation geom
               lmtri =  indmtr( i+1  ,j,k)
               if(neq_k.ne.3) then

                 do ne =1, neq_ij
                 ti( lmtr, ne)=ti( lmtri,ne)-  tj(lmtr,ne)+ tj(lmtrj,ne)
                 ti0(lmtr, ne)=ti0(lmtri,ne)- tj0(lmtr,ne)+tj0(lmtrj,ne)
                 enddo

               else
                 do ne =1, neq_ij
                 ti( lmtr, ne)= ti(lmtri,ne)- tj(lmtr,ne)+ tj(lmtrj,ne)
     &                                      - tk(lmtr,ne)+ tk(lmtrk,ne)
                 ti0(lmtr, ne)=ti0(lmtri,ne)-tj0(lmtr,ne)+ tj0(lmtrj,ne)
     &                                      -tk0(lmtr,ne)+ tk0(lmtrk,ne)
                 enddo
              endif

               vol(lmtr) = vol(lmtr0)

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
         do i= ind_dm(2)+1 , ind_dm(2)+nijk_mtr(4)

             jmin = max( ind_dm(3)-nijk_mtr(4) , j-1 )
             kmin = max( ind_dm(5)-nijk_mtr(5) , k-1 )
             jmax = min (ind_dm(4)+nijk_mtr(4) , j+1 )
             kmax = min (ind_dm(6)+nijk_mtr(5) , k+1 )

             ind0 = indcg(i, j , k)

             if(   degener(ind0).eq.1 )  then

               lmtr  =  indmtr( i                 , j, k)
               lmtr0 =  indmtr( 2*ind_dm(2)+ 1 -i , j, k)

               ind1 = indcg(i, jmin ,k  )
               if(degener(ind1).eq.1) then
                 do ne =1, neq_ij
                   tj( lmtr , ne) = tj( lmtr0 , ne)
                   tj0(lmtr , ne) = tj0(lmtr0 , ne)
                 enddo
               endif

               ind1 = indcg(i, j    ,kmin)
               if(degener(ind1).eq.1) then
                 do ne =1, neq_k
                   tk (lmtr  , ne)  = tk ( lmtr0 , ne)
                   tk0(lmtr  , ne)  = tk0( lmtr0, ne)
                 enddo
               endif

               lmtrj =  indmtr( i                 , j+1, k)
               lmtr0j=  indmtr( 2*ind_dm(2)+ 1 -i , j+1, k)

               ind1 = indcg(i,jmax,k)
               if(degener(ind1).eq.1) then

                 do ne =1, neq_ij
                   tj( lmtrj, ne) = tj( lmtr0j, ne)
                   tj0(lmtrj, ne) = tj0(lmtr0j, ne)
                 enddo
               endif

               lmtrk =  indmtr( i                 , j, k+1)
               lmtr0k=  indmtr( 2*ind_dm(2)+ 1 -i , j, k+1)

               ind1 = indcg(i,j,kmax)
               if(degener(ind1).eq.1) then

                 do ne =1, neq_k
                   tk(lmtrk , ne)  = tk( lmtr0k, ne)
                   tk0(lmtrk, ne)  = tk0(lmtr0k, ne)
                 enddo
               endif

               !!Reconstruction ti pour assurer conservation geom
               lmtri  =  indmtr( i+1              , j  , k)
               if(neq_k.ne.3) then

                 do ne =1, neq_ij
                 ti( lmtri, ne)=ti( lmtr,ne)+  tj(lmtr,ne)- tj(lmtrj,ne)
                 ti0(lmtri, ne)=ti0(lmtr,ne)+ tj0(lmtr,ne)-tj0(lmtrj,ne)
                 enddo

               else
                 do ne =1, neq_ij
                 ti( lmtri, ne)= ti(lmtr,ne)+ tj(lmtr,ne)- tj(lmtrj,ne)
     &                                      + tk(lmtr,ne)- tk(lmtrk,ne)
                 ti0(lmtri, ne)=ti0(lmtr,ne)+tj0(lmtr,ne)- tj0(lmtrj,ne)
     &                                      +tk0(lmtr,ne)- tk0(lmtrk,ne)
                 enddo
              endif

              vol(lmtr) = vol(lmtr0)

             endif

         enddo
         enddo
         enddo

      enddo !npass

      end

