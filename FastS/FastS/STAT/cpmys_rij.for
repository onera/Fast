c***********************************************************************
c     $Date: 2011-10-10 16:10:53 +0200 (lun 10 oct 2011) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cpmys_rij(ndom, ndimdx, ndimdx_mtr, ndimdmy,
     &                     neq, neq_my,neq_grad,neq_ij, neq_k,type_zone,
     &                     lthermique, lreynolds,
     &                     nijk_mtr, nijk, nijk_my,  ijkv,
     &                     moy_param,
     &                     ind_loop, 
     &                     gamma, cv, prandtl, 
     &                     rop,  xmut, 
     &                     ti, tj, tk, vol, ti_df, tj_df, tk_df,vol_df,
     &                     moy)
c***********************************************************************
c_U   USER : C. laurent
c
c     ACT
c_A    calcul grandeur moyenne pour bilan rij
c_A    commande de calcul des statistiques 
c
c     VAL
c
c     INP
c***********************************************************************
      implicit none

      logical lorder4
      INTEGER_E ndom,neq,neq_grad,neq_my,neq_ij, neq_k,
     & lthermique,lreynolds,type_zone,
     & ndimdx, ndimdmy,ndimdx_mtr, moy_param(*), ind_loop(6),
     & nijk_mtr(5), nijk(5), nijk_my(5), ijkv(3)
c
      REAL_E  rop(ndimdx, neq), xmut(ndimdx)
      REAL_E  moy(ndimdmy, neq_my), gamma, cv, prandtl

      REAL_E ti(ndimdx_mtr,neq_ij)   , tj(ndimdx_mtr,neq_ij),
     &       tk(ndimdx_mtr,neq_k)    , vol(ndimdx_mtr) 
      REAL_E ti_df(ndimdx_mtr,neq_ij), tj_df(ndimdx_mtr,neq_ij),
     &       tk_df(ndimdx_mtr,neq_k) , vol_df(ndimdx_mtr)

c Var loc
      INTEGER_E i,j,k,l,m,ideb,ifin,jdeb,jfin,kdeb,kfin,neq_tensrey,
     &        ind, ne,nd
      INTEGER_E  c1,c2,c3, eq0
      REAL_E    u1,u2,u3,p,mu,dukdxk, s11,s22,s33,s12,s23,s13
      REAL_E    du1dx1,du1dx2,du1dx3,du2dx1,du2dx2,du2dx3,cnm, cn
      REAL_E    du3dx1,du3dx2,du3dx3,cmu,nu,c4,cvgam,rho,rho_1,t,rg
      REAL_E    rho_cn,rou_cn,rov_cn,row_cn,p_cn,mu_cn,t_cn,nu_cn

#include "FastS/formule.h"
#include "FastS/formule_mtr.h"

      INTEGER_E indmy
      indmy(i_3,j_3,k_3) =  1
     &                   + (i_3+nijk_my(4)-1)*nijk_my(1)
     &                   + (j_3+nijk_my(4)-1)*nijk_my(2)
     &                   + (k_3+nijk_my(5)-1)*nijk_my(3)

      ! coeff pour moyenne glissante temporelle
      c1           = moy_param(3)                    !nbr echantillon instant N
      c2           = moy_param(3)-1                  !nbr echantillon instant N-1
      c3           = moy_param(6)                    !nbr de cellule homogene sur lequel on somme les variables dans cpmys
      cnm          = c2/float(c1)
      cn           = 1./float(c1*c3)

      rg       = cv*(gamma-1.)
      c4       = prandtl/cv/gamma

      neq_tensrey= 0
      if(lreynolds.eq.1) neq_tensrey = 80
 
      ideb= ind_loop(1)            
      jdeb= ind_loop(3)
      kdeb= ind_loop(5)
      ifin= ind_loop(2)
      jfin= ind_loop(4)
      kfin= ind_loop(6)

      if    ((moy_param(1)==1).and.(moy_param(2).eq.3)) then
        kdeb     = 1
        kfin     = ijkv(3)
      elseif((moy_param(1)==1).and.(moy_param(2).eq.2)) then
        jdeb     = 1
        jfin     = ijkv(2)
      elseif((moy_param(1)==1).and.(moy_param(2).eq.1)) then
        ideb     = 1
        ifin     = ijkv(1)
      elseif((moy_param(1)==2).and.(moy_param(2).eq.1)) then
        jdeb     = 1
        jfin     = ijkv(2)
        kdeb     = 1
        kfin     = ijkv(3)
      elseif((moy_param(1)==2).and.(moy_param(2).eq.2)) then
        ideb     = 1
        ifin     = ijkv(1)
        kdeb     = 1
        kfin     = ijkv(3)
      elseif((moy_param(1)==2).and.(moy_param(2).eq.3)) then
        ideb     = 1
        ifin     = ijkv(1)
        jdeb     = 1
        jfin     = ijkv(2)
      endif

c      write(*,'(a,6i5)')'loop_final',ideb,ifin,jdeb,jfin,kdeb,kfin
c      write(*,'(a,5i5)')'nijk',nijk
c      write(*,'(a,i5)')'ndimdmy',ndimdmy
c      write(*,'(a,5i5)')'nijk_my',nijk_my
c      write(*,'(a,3i5)')'neq_ij',neq_ij, neq_k
c      write(*,'(a,3i5)')'c1,c2',c1,c2,c3
c      write(*,'(a,2f18.9)')'cn,cm',cn*c3,cnm


!   Initialisation par la moyenne des n-1 echantillons
#ifndef E_PERMUT  
      do 10 ne=1,neq_my
        do 10 k= ind_loop(5), ind_loop(6)
        do 10 j= ind_loop(3), ind_loop(4)
!DEC$ IVDEP
        do 10 i= ind_loop(1), ind_loop(2)
         l =indmy( i, j, k)
         moy(l, ne) =   cnm*moy(l, ne)
 10   continue 
#else
      do 10 k= ind_loop(5), ind_loop(6)
      do 10 j= ind_loop(3), ind_loop(4)
      do 10 i= ind_loop(1), ind_loop(2)
         l =indmy( i, j, k)
!DEC$ IVDEP
        do 10 ne=1,neq_my
         moy(l, ne) =   cnm*moy(l, ne)
 10   continue 
#endif


      IF(lthermique.eq.0.and.lreynolds.eq.0) THEN  !moyenne favre basique: u et u^2 

         if(type_zone.ne.3) then  !dom 3D

             do k=kdeb,kfin
             do j=jdeb,jfin
!DEC$ IVDEP
             do i=ideb,ifin
#include      "FastS/STAT/cpmys_rij_ufavre_3d.for"
             end do
             end do
             end do

         else !dom 2d

             do j=jdeb,jfin
!DEC$ IVDEP
             do i=ideb,ifin
#include      "FastS/STAT/cpmys_rij_ufavre_2d.for"
             end do
             end do
         endif       !2d/3D 


      ELSEIF(lthermique.eq.0.and.lreynolds.eq.1) THEN  !moyenne u + bilan Eq transport rij

c         if(type_zone.ne.3) then  !dom 3D
c
c             do k=kdeb,kfin
c             do j=jdeb,jfin
c             do i=ideb,ifin
c#include      "FastS/STAT/cpmys_rij_ufavre_3d.for"
c#include      "FastS/STAT/cpmys_rij_tensrey_3d.for"
c             end do
c             end do
c             end do
c
c         else !dom 2d
c
c             do j=jdeb,jfin
c             do i=ideb,ifin
c#include      "FastS/STAT/cpmys_rij_ufavre_2d.for"
c#include      "FastS/STAT/cpmys_rij_tensrey_2d.for"
cc             end do
c             end do
c         endif       !2d/3D 

      ELSEIF(lthermique.eq.1.and.lreynolds.eq.0) THEN  !moyenne u + T' + roui'T' + (grad t')^2


c         if(type_zone.ne.3) then  !dom 3D
c
c             do k=kdeb,kfin
c             do j=jdeb,jfin
c             do i=ideb,ifin
c#include      "FastS/STAT/cpmys_rij_ufavre_3d.for"
c#include      "FastS/STAT/cpmys_rij_thermique_3d.for"
c             end do
c             end do
c             end do
c
c         else !dom 2d
c
c             do j=jdeb,jfin
c             do i=ideb,ifin
c#include      "FastS/STAT/cpmys_rij_ufavre_2d.for"
c#include      "FastS/STAT/cpmys_rij_thermique_2d.for"
c             end do
c             end do
c         endif       !2d/3D 
c
c
      ELSE       ! u + thermique + reynolds
c
c         if(type_zone.ne.3) then  !dom 3D
c
c             do k=kdeb,kfin
c             do j=jdeb,jfin
c             do i=ideb,ifin
c#include      "FastS/STAT/cpmys_rij_ufavre_3d.for"
c#include      "FastS/STAT/cpmys_rij_thermique_3d.for"
c#include      "FastS/STAT/cpmys_rij_tensrey_3d.for"
c             end do
c             end do
c             end do
c
c         else !dom 2d
c
c             do j=jdeb,jfin
c             do i=ideb,ifin
c#include      "FastS/STAT/cpmys_rij_ufavre_2d.for"
c#include      "FastS/STAT/cpmys_rij_thermique_2d.for"
c#include      "FastS/STAT/cpmys_rij_tensrey_2d.for"
c             end do
c             end do
c         endif       !2d/3D 
c
      ENDIF ! type moyenne
   
      end
