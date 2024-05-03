c***********************************************************************
c     $Date: 2010-07-12 18:57:35 +0200 (Mon, 12 Jul 2010) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu_luSA(ndom,neq,neq_ij,neq_k,neq_coe,
     &                   ndimdx, ndimdx_mtr, nijk, nijk_mtr,
     &                   imtr, gamma, cv,
     &                   ind_loop,
     &                   llower,
     &                   drodm,rop,
     &                   ti,tj,tk,
     &                   coe)
c***********************************************************************
c                              O N E R A
c
c_D   DATE_C/M : 1996
c
c_U   USER : DARRACQ 
c
c     ACT
c_A    Construction et inversion d'une matrice inferieure tetradiagonale par blocs
c      par une methode vectorisable.
c     VAL
c_V    Steady
c_V    Formulation LCI+Jameson-Turkel
c
c     INP
c_I    ndom,neq,drodm,coe
c
c     OUT
c     I/O
c_/    drodm,drodm
c***********************************************************************
      implicit none

      logical llower
      INTEGER_E ndom,neq,neq_ij,neq_k,neq_coe,ndimdx,ndimdx_mtr, imtr,
     & nijk(5), nijk_mtr(5), ind_loop(6)
 
      REAL_E gamma, cv
      REAL_E drodm(ndimdx,neq),coe(ndimdx,neq_coe),rop(ndimdx,neq)
      REAL_E ti(ndimdx_mtr,neq_ij),tj(ndimdx_mtr,neq_ij),
     &       tk(ndimdx_mtr,neq_k)

c Var loc
      INTEGER_E  inci,incj,inck,l,i,j,k,kdmax,kd,lmax,ll,ndo,
     & kddeb,kdfin,ipas,kfin,kdeb,jfin,jdeb,ifin,ideb,
     & l1,l2,lt,lt1,lt2,
     & inci2_mtr,incj2_mtr,inck2_mtr,inci_mtr,incj_mtr,inck_mtr

      REAL_E gam2,gam1,gamm1,cp,xal,diag,
     & b11,b12,b13,b14,b15,b21,b22,b23,b24,b25,b31,b32,b33,b34,b35,b41,
     & b42,b43,b44,b45,b51,b52,b53,b54,b55,
     & b1,b2,b3,b4,b5,b6,c1,c2,c3,c4,c5,c6,e1,e2,e3,e4,e5,e6,
     & signe,r,u,v,w,t,q2,h,ph2,qn,tcx,tcy,tcz

#include "FastS/formule.h"
#include "FastS/formule_mtr.h"

      gam1    = gamma
      gamm1   = gam1 - 1.
      cp      = gamma*cv
      gam2    = gamm1-1.
 
      if(llower) then
       inci      = 1
       incj      = nijk(1)
       inck      = nijk(1)*nijk(2)
       inci_mtr  = nijk_mtr(1)
       incj_mtr  = nijk_mtr(2)
       inck_mtr  = nijk_mtr(3)
       inci2_mtr = 0
       incj2_mtr = 0
       inck2_mtr = 0
       ipas      = 1
       signe     = 0.5
      else
       inci      =-1
       incj      =-nijk(1)
       inck      =-nijk(1)*nijk(2)
       inci_mtr  =-nijk_mtr(1)
       incj_mtr  =-nijk_mtr(2)
       inck_mtr  =-nijk_mtr(3)
       inci2_mtr = 2*inci_mtr
       incj2_mtr = 2*incj_mtr
       inck2_mtr = 2*inck_mtr
       ipas      =-1
       signe     =-0.5
      endif

      if(llower) then
        kdeb  = ind_loop(5)
        jdeb  = ind_loop(3)
        ideb  = ind_loop(1)
        kfin  = ind_loop(6)
        jfin  = ind_loop(4)
        ifin  = ind_loop(2)
      else
        kfin  = ind_loop(5)
        jfin  = ind_loop(3)
        ifin  = ind_loop(1)
        kdeb  = ind_loop(6)
        jdeb  = ind_loop(4)
        ideb  = ind_loop(2)
      endif
 
      IF(imtr.eq.0) THEN !domaine 3d general


      !!! coin (ideb,jdeb,kdb)

        l = inddm(ideb,jdeb,kdeb)
        lt= indmtr(ideb,jdeb,kdeb)

        diag      = 1./coe(l,5)

        drodm(l,1)= drodm(l,1)*diag
        drodm(l,2)= drodm(l,2)*diag
        drodm(l,3)= drodm(l,3)*diag
        drodm(l,4)= drodm(l,4)*diag
        drodm(l,5)= drodm(l,5)*diag
        drodm(l,6)= drodm(l,6)/coe(l,6)

      !!! ligne (jdeb,kdeb) dans le plan kdeb

        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          lt     = indmtr(i,jdeb,kdeb)

#include "FastS/Compute/LU/lu_i_3dfull.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ b1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ b2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ b3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ b4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ b5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ b6*xal)/coe(l,6)
        enddo

      !!! ligne (ideb,kdeb) dans le plan kdeb

        do j= jdeb+ipas,jfin,ipas

          l      = inddm(ideb,j,kdeb)
          lt     = indmtr(ideb,j,kdeb)

#include "FastS/Compute/LU/lu_j_3dfull.for"
          c6  = (qn+coe(l1,3)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ c1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ c2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ c3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ c4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ c5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ c6*xal)/coe(l,6)
        enddo

        !!! Fin plan kdeb 

        do j= jdeb+ipas,jfin,ipas
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,j,kdeb)
          lt     = indmtr(i,j,kdeb)

#include "FastS/Compute/LU/lu_i_3dfull.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_3dfull.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+c1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+c2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+c3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(b4+c4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+c5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+c6)*xal)/coe(l,6)
        enddo
        enddo

        !!! ligne (ideb,jdb)  dans le plan jdeb

        do k= kdeb+ipas,kfin,ipas

          l      = inddm(ideb,jdeb,k)
          lt     = indmtr(ideb,jdeb,k)

#include "FastS/Compute/LU/lu_k_3dfull.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ e1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ e2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ e3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ e4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ e5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ e6*xal)/coe(l,6)
        enddo

        !!! Fin plan jdeb 
        do k= kdeb+ipas,kfin,ipas
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,k)
          lt     = indmtr(i,jdeb,k)

#include "FastS/Compute/LU/lu_i_3dfull.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dfull.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal     = coe(l,1)*signe*0.5
          diag    = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+e1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+e2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+e3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(b4+e4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+e5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+e6)*xal)/coe(l,6)
        enddo
        enddo

        !!! Fin plan ideb 
        do k= kdeb+ipas,kfin,ipas
        do j= jdeb+ipas,jfin,ipas

          l      = inddm(ideb,j,k)
          lt     = indmtr(ideb,j,k)

#include "FastS/Compute/LU/lu_j_3dfull.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dfull.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(c1+e1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(c2+e2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(c3+e3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(c4+e4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(c5+e5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(c6+e6)*xal)/coe(l,6)
        enddo
        enddo

        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas
        do  j= jdeb+ipas,jfin,ipas
        do  i= ideb+ipas,ifin,ipas
      
          l  = inddm(i,j,k)
          lt = indmtr(i,j,k)

#include "FastS/Compute/LU/lu_i_3dfull.for"
         b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_3dfull.for"
         c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dfull.for"
         e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

         xal    = coe(l,1)*signe*0.5
         diag   = 1./coe(l,5)
c
         drodm(l,1)= (drodm(l,1)+(b1+c1+e1)*xal)*diag
         drodm(l,2)= (drodm(l,2)+(b2+c2+e2)*xal)*diag
         drodm(l,3)= (drodm(l,3)+(b3+c3+e3)*xal)*diag
         drodm(l,4)= (drodm(l,4)+(b4+c4+e4)*xal)*diag
         drodm(l,5)= (drodm(l,5)+(b5+c5+e5)*xal)*diag
         drodm(l,6)= (drodm(l,6)+(b6+c6+e6)*xal)/coe(l,6) 
        enddo
        enddo
        enddo

      ELSEIF(imtr.eq.1) THEN !maillage 3d k homogene:

      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)

        l = inddm(ideb,jdeb,kdeb)
        lt= indmtr(ideb,jdeb,kdeb)

        diag      = 1./coe(l,5)

        drodm(l,1)= drodm(l,1)*diag
        drodm(l,2)= drodm(l,2)*diag
        drodm(l,3)= drodm(l,3)*diag
        drodm(l,4)= drodm(l,4)*diag
        drodm(l,5)= drodm(l,5)*diag
        drodm(l,6)= drodm(l,6)/coe(l,6)

      !!! ligne (jdeb,kdeb) dans le plan kdeb

        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          lt     = indmtr(i,jdeb,kdeb)

#include "FastS/Compute/LU/lu_i_3dhomogene.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ b1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ b2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ b3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ b4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ b5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ b6*xal)/coe(l,6)
        enddo

      !!! ligne (ideb,kdeb) dans le plan kdeb

        do j= jdeb+ipas,jfin,ipas

          l      = inddm(ideb,j,kdeb)
          lt     = indmtr(ideb,j,kdeb)

#include "FastS/Compute/LU/lu_j_3dhomogene.for"
          c6  = (qn+coe(l1,3)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ c1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ c2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ c3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ c4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ c5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ c6*xal)/coe(l,6)
        enddo

        !!! Fin plan kdeb 

        do j= jdeb+ipas,jfin,ipas
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,j,kdeb)
          lt     = indmtr(i,j,kdeb)

#include "FastS/Compute/LU/lu_i_3dhomogene.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_3dhomogene.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+c1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+c2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+c3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(b4+c4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+c5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+c6)*xal)/coe(l,6)
        enddo
        enddo

        !!! ligne (ideb,jdb)  dans le plan jdeb

        do k= kdeb+ipas,kfin,ipas

          l      = inddm(ideb,jdeb,k)
          lt     = indmtr(ideb,jdeb,k)

#include "FastS/Compute/LU/lu_k_3dhomogene.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ e1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ e2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ e3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ e4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ e5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ e6*xal)/coe(l,6)
        enddo

        !!! Fin plan jdeb 
        do k= kdeb+ipas,kfin,ipas
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,k)
          lt     = indmtr(i,jdeb,k)

#include "FastS/Compute/LU/lu_i_3dhomogene.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dhomogene.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal     = coe(l,1)*signe*0.5
          diag    = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+e1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+e2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+e3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(b4+e4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+e5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+e6)*xal)/coe(l,6)
        enddo
        enddo

        !!! Fin plan ideb 
        do k= kdeb+ipas,kfin,ipas
        do j= jdeb+ipas,jfin,ipas

          l      = inddm(ideb,j,k)
          lt     = indmtr(ideb,j,k)

#include "FastS/Compute/LU/lu_j_3dhomogene.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dhomogene.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(c1+e1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(c2+e2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(c3+e3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(c4+e4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(c5+e5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(c6+e6)*xal)/coe(l,6)
        enddo
        enddo

        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas
        do  j= jdeb+ipas,jfin,ipas
        do  i= ideb+ipas,ifin,ipas
      
         l = inddm(i,j,k)
         lt= indmtr(i,j,k)

#include "FastS/Compute/LU/lu_i_3dhomogene.for"
         b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_3dhomogene.for"
         c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dhomogene.for"
         e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

         xal    = coe(l,1)*signe*0.5
         diag   = 1./coe(l,5)
c
         drodm(l,1)= (drodm(l,1)+(b1+c1+e1)*xal)*diag
         drodm(l,2)= (drodm(l,2)+(b2+c2+e2)*xal)*diag
         drodm(l,3)= (drodm(l,3)+(b3+c3+e3)*xal)*diag
         drodm(l,4)= (drodm(l,4)+(b4+c4+e4)*xal)*diag
         drodm(l,5)= (drodm(l,5)+(b5+c5+e5)*xal)*diag
         drodm(l,6)= (drodm(l,6)+(b6+c6+e6)*xal)/coe(l,6) 
        enddo
        enddo
        enddo

      ELSEIF(imtr.eq.2) THEN !maillage 3d cartesien


       tcx = ti(1,1)
       tcy = tj(1,1)
       tcz = tk(1,1)
      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)

        l = inddm(ideb,jdeb,kdeb)
        lt= indmtr(ideb,jdeb,kdeb)

        diag      = 1./coe(l,5)

        drodm(l,1)= drodm(l,1)*diag
        drodm(l,2)= drodm(l,2)*diag
        drodm(l,3)= drodm(l,3)*diag
        drodm(l,4)= drodm(l,4)*diag
        drodm(l,5)= drodm(l,5)*diag
        drodm(l,6)= drodm(l,6)/coe(l,6)

      !!! ligne (jdeb,kdeb) dans le plan kdeb

        do i= ideb+ipas,ifin,ipas

          l     = inddm(i,jdeb,kdeb)
          lt    = indmtr(i,jdeb,kdeb)

#include "FastS/Compute/LU/lu_i_3dcart.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ b1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ b2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ b3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ b4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ b5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ b6*xal)/coe(l,6)
        enddo

      !!! ligne (ideb,kdeb) dans le plan kdeb

        do j= jdeb+ipas,jfin,ipas

          l     = inddm(ideb,j,kdeb)
          lt    = indmtr(ideb,j,kdeb)

#include "FastS/Compute/LU/lu_j_3dcart.for"
          c6  = (qn+coe(l1,3)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ c1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ c2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ c3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ c4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ c5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ c6*xal)/coe(l,6)
        enddo

        !!! Fin plan kdeb 

        do j= jdeb+ipas,jfin,ipas
        do i= ideb+ipas,ifin,ipas

          l     = inddm(i,j,kdeb)
          lt    = indmtr(i,j,kdeb)

#include "FastS/Compute/LU/lu_i_3dcart.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_3dcart.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+c1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+c2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+c3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(b4+c4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+c5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+c6)*xal)/coe(l,6)
        enddo
        enddo

        !!! ligne (ideb,jdb)  dans le plan jdeb

        do k= kdeb+ipas,kfin,ipas

          l     = inddm(ideb,jdeb,k)
          lt    = indmtr(ideb,jdeb,k)

#include "FastS/Compute/LU/lu_k_3dcart.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ e1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ e2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ e3*xal)*diag
          drodm(l,4)= (drodm(l,4)+ e4*xal)*diag
          drodm(l,5)= (drodm(l,5)+ e5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ e6*xal)/coe(l,6)
        enddo

        !!! Fin plan jdeb 
        do k= kdeb+ipas,kfin,ipas
        do i= ideb+ipas,ifin,ipas

          l     = inddm(i,jdeb,k)
          lt    = indmtr(i,jdeb,k)

#include "FastS/Compute/LU/lu_i_3dcart.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dcart.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal     = coe(l,1)*signe*0.5
          diag    = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+e1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+e2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+e3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(b4+e4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+e5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+e6)*xal)/coe(l,6)
        enddo
        enddo

        !!! Fin plan ideb 
        do k= kdeb+ipas,kfin,ipas
        do j= jdeb+ipas,jfin,ipas

          l     = inddm(ideb,j,k)
          lt    = indmtr(ideb,j,k)

#include "FastS/Compute/LU/lu_j_3dcart.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dcart.for"
          e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

          xal    = coe(l,1)*signe*0.5
          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(c1+e1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(c2+e2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(c3+e3)*xal)*diag
          drodm(l,4)= (drodm(l,4)+(c4+e4)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(c5+e5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(c6+e6)*xal)/coe(l,6)
        enddo
        enddo

        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas
        do  j= jdeb+ipas,jfin,ipas
        do  i= ideb+ipas,ifin,ipas
      
         l  = inddm(i,j,k)
         lt = indmtr(i,j,k)

#include "FastS/Compute/LU/lu_i_3dcart.for"
         b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_3dcart.for"
         c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_k_3dcart.for"
         e6 = (qn+coe(l1,4)*signe)*drodm(l1,6)

         xal    = coe(l,1)*signe*0.5
         diag   = 1./coe(l,5)
c
         drodm(l,1)= (drodm(l,1)+(b1+c1+e1)*xal)*diag
         drodm(l,2)= (drodm(l,2)+(b2+c2+e2)*xal)*diag
         drodm(l,3)= (drodm(l,3)+(b3+c3+e3)*xal)*diag
         drodm(l,4)= (drodm(l,4)+(b4+c4+e4)*xal)*diag
         drodm(l,5)= (drodm(l,5)+(b5+c5+e5)*xal)*diag
         drodm(l,6)= (drodm(l,6)+(b6+c6+e6)*xal)/coe(l,6) 
        enddo
        enddo
        enddo


      ELSE !2d

      !!! coin (ideb,jdeb,kdb)

        l = inddm(ideb,jdeb,1)
        lt= indmtr(ideb,jdeb,1)

        diag      = 1./coe(l,5)

        drodm(l,1)= drodm(l,1)*diag
        drodm(l,2)= drodm(l,2)*diag
        drodm(l,3)= drodm(l,3)*diag
        drodm(l,4)= drodm(l,4)*diag
        drodm(l,5)= drodm(l,5)*diag
        drodm(l,6)= drodm(l,6)/coe(l,6)

      !!! ligne (jdeb,kdeb) dans le plan kdeb

        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,1)
          lt     = indmtr(i,jdeb,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_2d.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)

          diag    = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ b1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ b2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ b3*xal)*diag
          drodm(l,5)= (drodm(l,5)+ b5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ b6*xal)/coe(l,6)
        enddo

      !!! ligne (ideb,kdeb) dans le plan kdeb

      do j= jdeb+ipas,jfin,ipas

          l      = inddm(ideb,j,1)
          lt     = indmtr(ideb,j,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_2d.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)

          diag    = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+ c1*xal)*diag
          drodm(l,2)= (drodm(l,2)+ c2*xal)*diag
          drodm(l,3)= (drodm(l,3)+ c3*xal)*diag
          drodm(l,5)= (drodm(l,5)+ c5*xal)*diag
          drodm(l,6)= (drodm(l,6)+ c6*xal)/coe(l,6)
        enddo

      !!! Fin plan kdeb

      do j= jdeb+ipas,jfin,ipas
      do i= ideb+ipas,ifin,ipas

          l      = inddm(i,j,1)
          lt     = indmtr(i,j,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_2d.for"
          b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)
#include "FastS/Compute/LU/lu_j_2d.for"
          c6 = (qn+coe(l1,3)*signe)*drodm(l1,6)

          diag   = 1./coe(l,5)

          drodm(l,1)= (drodm(l,1)+(b1+c1)*xal)*diag
          drodm(l,2)= (drodm(l,2)+(b2+c2)*xal)*diag
          drodm(l,3)= (drodm(l,3)+(b3+c3)*xal)*diag
          drodm(l,5)= (drodm(l,5)+(b5+c5)*xal)*diag
          drodm(l,6)= (drodm(l,6)+(b6+c6)*xal)/coe(l,6)
        enddo
        enddo

      ENDIF
 
      end
