c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine skmtr_df(nijk_xyz, nijk_mtr,
     &                    ndimt_xyz, x, y, z,
     &                    tmpi,tmpj,tmpk,tmpi2,tmpj2,tmpk2, mtr,
     &                    neq_ij, neq_k,
     &                    ndimdx_mtr, ti,tj,tk, detjbar)
c***********************************************************************
c     ACT
c_A    Calcul de la metrique sur l'ensemble d'une configuration
c      dans le cas d'un schema aux differences finies.
c      M. Terracol,   Novembre 2001
c
c     VAL
c_V    Steady/Unsteady !
c
c     INP
c_I    ndom     : numero du domaine courant
c_I    x,y,z    : coordonnees des centres de cellules
c
c     OUT
c_O    t(ijk)d  : composantes de J/det(J)
c_O    detjbar  : det(Jbar) = 1/det(J)
c***********************************************************************
      implicit none

      INTEGER_E nijk_mtr(5),nijk_xyz(5),ndimt_xyz,ndimdx_mtr,ndimf,
     & neq_ij, neq_k, ind_ordre

      REAL_E x(ndimt_xyz),y(ndimt_xyz),z(ndimt_xyz)
      REAL_E tmpi(ndimt_xyz) ,tmpj(ndimt_xyz) ,tmpk(ndimt_xyz)
      REAL_E tmpi2(ndimt_xyz),tmpj2(ndimt_xyz),tmpk2(ndimt_xyz)
      REAL_E mtr(ndimt_xyz,9)

      REAL_E ti(ndimdx_mtr,neq_ij),tj(ndimdx_mtr,neq_ij),
     &       tk(ndimdx_mtr,neq_k), detjbar(ndimdx_mtr)

C Var loc
      INTEGER_E l1,l2,l3,l4,l5,l6,l7,l8,kficc,iormtr,iormtr2,
     & i1,i2,j1,j2,
     & i,ifin,j,jfin,k,kfin,k1,k2,deb,ideb,jdeb,kdeb,idir,l,ificc,
     & ific_min,incr,ijk_bord,ishift,ne, iordre_mtr,iorder_sch,iv,jv,kv
      REAL_E dz,jbar11,jbar12,jbar13,jbar21,jbar22,jbar23,jbar31,
     & jbar32,jbar33,x1,x2,x3

#include "FastS/formule_mtr.h"

C    adresse point courant pour tableau x ou de la taille d'un domaine 
      INTEGER_E indcg, i_1,j_1,k_1

      indcg(i_1,j_1,k_1) =  1
     &                   + (i_1+nijk_xyz(4)-1)
     &                   + (j_1+nijk_xyz(4)-1)*nijk_xyz(1)
     &                   + (k_1+nijk_xyz(5)-1)*nijk_xyz(1)*nijk_xyz(2)

      iordre_mtr = 6
      iorder_sch = 6
      ind_ordre  = 0

      iv = nijk_xyz(1)-2*nijk_xyz(4) -1
      jv = nijk_xyz(2)-2*nijk_xyz(4) -1
      kv = nijk_xyz(3)-2*nijk_xyz(5) -1

      if(nijk_xyz(5).eq.0) kv = 1

      !write(*,'(a,5i5)')'nijk_xyz',nijk_xyz
c.... flag associe a l'option DF 2D
 
      if(ind_ordre.eq.1) then
        kficc = min(nijk_mtr(4),nijk_mtr(5))
        iormtr  = 2
        iormtr2 = 2
        i1 = 1  - nijk_mtr(4)
        i2 = iv + nijk_mtr(4)
        j1 = 1  - nijk_mtr(4)
        j2 = jv + nijk_mtr(4)
        k1 = 1  - kficc
        k2 = kv + kficc
      elseif(ind_ordre.eq.2) then
        kficc = min(2,nijk_mtr(5))
        iormtr  = 4 
        iormtr2 = 2
        i1 =    - 1
        i2 = iv + 2
        j1 =    - 1
        j2 = jv + 2
        k1 = 1  - kficc
        k2 = kv + kficc
      elseif(ind_ordre.eq.3) then
        kficc   = min(1,nijk_mtr(5))
        iormtr  = 6 
        iormtr2 = 2
        i1 = 0
        i2 = iv + 1
        j1 = 0
        j2 = jv + 1
        k1 = 1  - kficc
        k2 = kv + kficc
      else
        kficc = min(nijk_mtr(4),nijk_mtr(5))
        iormtr  = iordre_mtr
        iormtr2 = iordre_mtr
        i1 = 1  - nijk_mtr(4)
        i2 = iv + nijk_mtr(4)
        j1 = 1  - nijk_mtr(4)
        j2 = jv + nijk_mtr(4)
        k1 = 1  - kficc
        k2 = kv + kficc
      endif  

      write(*,'(a,6i5)')'i1,..k2',i1,i2,j1,j2,k1,k2
       
      !call inicfdf( ndimt_xyz, ndimf, cfdf,ordreloc,iordre_mtr)

      !on calcul le pas d'espace en z pour domaine 2d
      dz = (z(indcg(1,1,2))-z(indcg(1,1,1))) 

      ific_min = iorder_sch/2+iormtr2/2
      if(ific_min.gt.nijk_xyz(4)) then
        call error('skmtr_df$',301,1)
      endif

      ! Calcul des coordonnee
      deb  = 1 -nijk_xyz(4)
      kdeb = 1 -nijk_xyz(5)
      ifin = iv+nijk_xyz(4)
      jfin = jv+nijk_xyz(4)
      kfin = kv+nijk_xyz(5)
 
      do k = kdeb,kfin
      do j = deb,jfin
      do i = deb,ifin
            l1 = indcg(i,j,k)
            l2 = indcg(i+1,j,k)
            l3 = indcg(i+1,j+1,k)
            l4 = indcg(i+1,j+1,k+1)
            l5 = indcg(i+1,j,k+1)
            l6 = indcg(i,j+1,k+1)
            l7 = indcg(i,j+1,k)
            l8 = indcg(i,j,k+1)
            l  = indcg(i,j,k)

            tmpi(l) = (x(l1)+x(l2)+x(l3)+x(l4)
     &                 +x(l5)+x(l6)+x(l7)+x(l8))/8. 
            tmpj(l) = (y(l1)+y(l2)+y(l3)+y(l4)
     &                 +y(l5)+y(l6)+y(l7)+y(l8))/8. 
            tmpk(l) = (z(l1)+z(l2)+z(l3)+z(l4)
     &                 +z(l5)+z(l6)+z(l7)+z(l8))/8. 
      enddo
      enddo
      enddo

      call derdf_xyz(nijk_xyz, ndimt_xyz, tmpi,tmpi2 , iormtr, 1)
      call derdf_xyz(nijk_xyz, ndimt_xyz, tmpi,tmpj2 , iormtr, 2)
      call derdf_xyz(nijk_xyz, ndimt_xyz, tmpi,tmpk2 , iormtr, 3)

      do l =1,ndimt_xyz
          mtr(l,1) = tmpi2(l)
          mtr(l,4) = tmpj2(l)
          mtr(l,7) = tmpk2(l)
      enddo

      call derdf_xyz(nijk_xyz, ndimt_xyz, tmpj,tmpi2 , iormtr, 1)
      call derdf_xyz(nijk_xyz, ndimt_xyz, tmpj,tmpj2 , iormtr, 2)
      call derdf_xyz(nijk_xyz, ndimt_xyz, tmpj,tmpk2 , iormtr, 3)

      do l =1,ndimt_xyz
          mtr(l,2) = tmpi2(l)
          mtr(l,5) = tmpj2(l)
          mtr(l,8) = tmpk2(l)
      enddo

      if (nijk_mtr(3).ne.1) then !dom3d
        call derdf_xyz(nijk_xyz, ndimt_xyz, tmpk,tmpi2 , iormtr, 1)
        call derdf_xyz(nijk_xyz, ndimt_xyz, tmpk,tmpj2 , iormtr, 2)
        call derdf_xyz(nijk_xyz, ndimt_xyz, tmpk,tmpk2 , iormtr, 3)
      else
          tmpi2(:) = 0.
          tmpj2(:) = 0.
          tmpk2(:) = dz
      endif

      do l =1,ndimt_xyz
          mtr(l,3) = tmpi2(l)
          mtr(l,6) = tmpj2(l)
          mtr(l,9) = tmpk2(l)
      enddo

      !Par la suite:
      !Jbar_ij    = d(X_j)/d(Xsi_i)
      !t(j)d(:,i) = J_ij/det(J)      = d(Xsi_j)/d(X_i) / det(J)
 
      do k = k1,k2
      do j = j1,j2
      do i = i1,i2

         l1      = indcg(i,j,k)
         l       = indmtr(i,j,k)

          jbar11 = mtr(l1,1)
          jbar12 = mtr(l1,2)
          jbar13 = mtr(l1,3)
          jbar21 = mtr(l1,4)
          jbar22 = mtr(l1,5)
          jbar23 = mtr(l1,6)
          jbar31 = mtr(l1,7)
          jbar32 = mtr(l1,8)
          jbar33 = mtr(l1,9)


          detjbar(l)  = max((
     &              jbar11*(jbar22*jbar33-jbar32*jbar23)
     &            - jbar12*(jbar21*jbar33-jbar31*jbar23)
     &            + jbar13*(jbar21*jbar32-jbar31*jbar22)),
     &            1.e-15)
         end do
         end do
         end do

        do l = 1,ndimt_xyz

          jbar11 = mtr(l,1)
          jbar12 = mtr(l,2)
          jbar13 = mtr(l,3)
          jbar21 = mtr(l,4)
          jbar22 = mtr(l,5)
          jbar23 = mtr(l,6)
          jbar31 = mtr(l,7)
          jbar32 = mtr(l,8)
          jbar33 = mtr(l,9)

          x1 = tmpi(l)
          x2 = tmpj(l)
          x3 = tmpk(l)

          tmpi(l)  = jbar22
          tmpk(l)  = jbar12

           mtr(l,1) = x1*jbar13
           mtr(l,2) = x1*jbar23
           mtr(l,3) = x1*jbar33
           mtr(l,4) = x2*jbar11
           mtr(l,5) = x2*jbar21
           mtr(l,6) = x2*jbar31
           mtr(l,7) = x3*jbar12
           mtr(l,8) = x3*jbar22
           mtr(l,9) = x3*jbar32
         end do

         !Terme 11
         if (nijk_mtr(3).ne.1) then !dom3d
           call derdf_xyz(nijk_xyz,ndimt_xyz, mtr(:,8),tmpi ,iormtr,3)
         else
           tmpi(:) = dz * tmpi(:)
         endif

         call derdf_xyz(nijk_xyz, ndimt_xyz, mtr(:,9),tmpj ,iormtr,2)
         do k = k1,k2
         do j = j1,j2
         do i = i1,i2
           l1      = indcg(i,j,k)
           l       = indmtr(i,j,k)
           ti(l,1) = tmpi(l1)-tmpj(l1)
         end do
         end do
         end do

         !Terme 12
         call derdf_xyz(nijk_xyz, ndimt_xyz, mtr(:,2),tmpi , iormtr, 3)
         call derdf_xyz(nijk_xyz, ndimt_xyz, mtr(:,3),tmpj , iormtr, 2)
         do k = k1,k2
         do j = j1,j2
         do i = i1,i2
           l1      = indcg(i,j,k)
           l       = indmtr(i,j,k)
           ti(l,2) = tmpi(l1)-tmpj(l1)
         end do
         end do
         end do

         !Terme 13
         if(nijk_mtr(3).ne.1.and.neq_k.eq.3) then
           call derdf_xyz(nijk_xyz,ndimt_xyz, mtr(:,5),tmpi ,iormtr, 3)
           call derdf_xyz(nijk_xyz,ndimt_xyz, mtr(:,6),tmpj ,iormtr, 2)
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
             l1      = indcg(i,j,k)
             l       = indmtr(i,j,k)
             ti(l,3) = tmpi(l1)-tmpj(l1)
           end do
           end do
           end do
         endif!2d/3d

         !Terme 21
         call derdf_xyz(nijk_xyz, ndimt_xyz, mtr(:,9),tmpi , iormtr, 1)
         if(nijk_mtr(3).ne.1) then
           call derdf_xyz(nijk_xyz,ndimt_xyz, mtr(:,7),tmpj ,iormtr, 3)
         else
           tmpj(:) = dz*tmpk(:)
         endif
         do k = k1,k2
         do j = j1,j2
         do i = i1,i2
             l1      = indcg(i,j,k)
             l       = indmtr(i,j,k)
             tj(l,1) = tmpi(l1)-tmpj(l1)
         end do
         end do
         end do

         !Terme 22
         call derdf_xyz(nijk_xyz, ndimt_xyz, mtr(:,3),tmpi , iormtr, 1)
         call derdf_xyz(nijk_xyz, ndimt_xyz, mtr(:,1),tmpj , iormtr, 3)
          do k = k1,k2
          do j = j1,j2
          do i = i1,i2
             l1      = indcg(i,j,k)
             l       = indmtr(i,j,k)
             tj(l,2) = tmpi(l1)-tmpj(l1)
          end do
          end do
          end do

         !Terme 23
         if(nijk_mtr(3).ne.1.and.neq_k.eq.3) then
           call derdf_xyz(nijk_xyz,ndimt_xyz, mtr(:,6),tmpi ,iormtr, 1)
           call derdf_xyz(nijk_xyz,ndimt_xyz, mtr(:,4),tmpj ,iormtr, 3)
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
             l1      = indcg(i,j,k)
             l       = indmtr(i,j,k)
             tj(l,3) = tmpi(l1)-tmpj(l1)
           end do
           end do
           end do
         endif !2d/3d

         if(nijk_mtr(3).ne.1) then
           if(neq_k.eq.3) then !3d full  
           !Terme 31
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,7),tmpi ,iormtr,2)
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,8),tmpj ,iormtr,1)
             do k = k1,k2
             do j = j1,j2
             do i = i1,i2
               l1      = indcg(i,j,k)
               l       = indmtr(i,j,k)
               tk(l,1) = tmpi(l1)-tmpj(l1)
             end do
             end do
             end do
             !Terme 32
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,1),tmpi ,iormtr,2)
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,2),tmpj ,iormtr,1)
             do k = k1,k2
             do j = j1,j2
             do i = i1,i2
               l1      = indcg(i,j,k)
               l       = indmtr(i,j,k)
               tk(l,2) = tmpi(l1)-tmpj(l1)
             end do
             end do
             end do
             !Terme 33
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,4),tmpi ,iormtr,2)
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,5),tmpj ,iormtr,1)
             do k = k1,k2
             do j = j1,j2
             do i = i1,i2
               l1      = indcg(i,j,k)
               l       = indmtr(i,j,k)
               tk(l,3) = tmpi(l1)-tmpj(l1)
             end do
             end do
             end do

           elseif(neq_k.eq.1) then !3d  homogene 

            ! Terme 33
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,4),tmpi ,iormtr,2)
             call derdf_xyz(nijk_xyz,ndimt_xyz,mtr(:,5),tmpj ,iormtr,1)
             do k = k1,k2
             do j = j1,j2
             do i = i1,i2
               l1      = indcg(i,j,k)
               l       = indmtr(i,j,k)
               tk(l,1) = tmpi(l1)-tmpj(l1)
             end do
             end do
             end do
            else
             continue
            endif !neq_k=1,2 ou 3
          endif !2d/3d

      end
