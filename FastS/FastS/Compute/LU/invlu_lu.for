c***********************************************************************
c     $Date: 2010-07-12 18:57:35 +0200 (Mon, 12 Jul 2010) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu_lu(ndom,neq,neq_ij,neq_k,neq_coe,
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
     & inci2,incj2,inck2,l1,l2,lt,lt1,lt2,
     & inci2_mtr,incj2_mtr,inck2_mtr,inci_mtr,incj_mtr,inck_mtr

      REAL_E gam2,gam1,gamm1,cp,xal,diag,
     & b11,b12,b13,b14,b15,b21,b22,b23,b24,b25,b31,b32,b33,b34,b35,b41,
     & b42,b43,b44,b45,b51,b52,b53,b54,b55,
     & b1,b2,b3,b4,b5,c1,c2,c3,c4,c5,e1,e2,e3,e4,e5,
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
       inci2     = 0
       incj2     = 0
       inck2     = 0
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
       inci2     =-2*inci 
       incj2     =-2*incj
       inck2     =-2*inck
       inci2_mtr =-2*inci_mtr
       incj2_mtr =-2*incj_mtr
       inck2_mtr =-2*inck_mtr
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

#include "FastS/Compute/LU/lu_dinv.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          lt     = indmtr(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dfull.for"
#include "FastS/Compute/LU/lu_dinv.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =  inddm(ideb,j,kdeb)
          lt     = indmtr(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dfull.for"
#include "FastS/Compute/LU/lu_dinv.for"

          do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,kdeb)
            lt     = indmtr(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dfull.for"
#include    "FastS/Compute/LU/lu_j_3dfull.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
        enddo
        !!! plan kdeb termine
        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          l      =  inddm(ideb,jdeb,k)
          lt     = indmtr(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dfull.for"
#include "FastS/Compute/LU/lu_dinv.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =  inddm(i,jdeb,k)
             lt     = indmtr(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dfull.for"
#include    "FastS/Compute/LU/lu_k_3dfull.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =  inddm(ideb,j,k)
             lt     = indmtr(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dfull.for"
#include    "FastS/Compute/LU/lu_k_3dfull.for"
#include    "FastS/Compute/LU/lu_dinv.for"

             do  i= ideb+ipas,ifin,ipas
      
               l = inddm(i,j,k)
               lt= indmtr(i,j,k)

                xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dfull.for"
#include       "FastS/Compute/LU/lu_j_3dfull.for"
#include       "FastS/Compute/LU/lu_k_3dfull.for"
#include       "FastS/Compute/LU/lu_dinv.for"
             enddo
          enddo
        enddo


      ELSEIF(imtr.eq.1) THEN !maillage 3d k homogene:

      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)
        l = inddm(ideb,jdeb,kdeb)

#include "FastS/Compute/LU/lu_dinv.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          lt     = indmtr(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dhomogene.for"
#include "FastS/Compute/LU/lu_dinv.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =  inddm(ideb,j,kdeb)
          lt     = indmtr(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dhomogene.for"
#include "FastS/Compute/LU/lu_dinv.for"

          do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,kdeb)
            lt     = indmtr(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dhomogene.for"
#include    "FastS/Compute/LU/lu_j_3dhomogene.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
        enddo
        !!! plan kdeb termine
        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          l      =  inddm(ideb,jdeb,k)
          lt     = indmtr(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dhomogene.for"
#include "FastS/Compute/LU/lu_dinv.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =  inddm(i,jdeb,k)
             lt     = indmtr(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dhomogene.for"
#include    "FastS/Compute/LU/lu_k_3dhomogene.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =  inddm(ideb,j,k)
             lt     = indmtr(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dhomogene.for"
#include    "FastS/Compute/LU/lu_k_3dhomogene.for"
#include    "FastS/Compute/LU/lu_dinv.for"

             do  i= ideb+ipas,ifin,ipas
      
               l = inddm(i,j,k)
               lt= indmtr(i,j,k)

                xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dhomogene.for"
#include       "FastS/Compute/LU/lu_j_3dhomogene.for"
#include       "FastS/Compute/LU/lu_k_3dhomogene.for"
#include       "FastS/Compute/LU/lu_dinv.for"
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

#include "FastS/Compute/LU/lu_dinv.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dcart.for"
#include "FastS/Compute/LU/lu_dinv.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =  inddm(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dcart.for"
#include "FastS/Compute/LU/lu_dinv.for"

          do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dcart.for"
#include    "FastS/Compute/LU/lu_j_3dcart.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
        enddo
        !!! plan kdeb termine
        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          l      =  inddm(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dcart.for"
#include "FastS/Compute/LU/lu_dinv.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =  inddm(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dcart.for"
#include    "FastS/Compute/LU/lu_k_3dcart.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =  inddm(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dcart.for"
#include    "FastS/Compute/LU/lu_k_3dcart.for"
#include    "FastS/Compute/LU/lu_dinv.for"

             do  i= ideb+ipas,ifin,ipas
      
               l = inddm(i,j,k)

                xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dcart.for"
#include       "FastS/Compute/LU/lu_j_3dcart.for"
#include       "FastS/Compute/LU/lu_k_3dcart.for"
#include       "FastS/Compute/LU/lu_dinv.for"
             enddo
          enddo
        enddo



      ELSE !2d

      !!! coin (ideb,jdeb,kdb)

        l = inddm(ideb,jdeb,1)
#include       "FastS/Compute/LU/lu_dinv_2d.for"

      !!! ligne (jdeb,kdeb) dans le plan kdeb

        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,1)
          lt     = indmtr(i,jdeb,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_2d.for"
#include "FastS/Compute/LU/lu_dinv_2d.for"
        enddo

      !!! ligne (ideb,kdeb) dans le plan kdeb

      do j= jdeb+ipas,jfin,ipas

          l      = inddm(ideb,j,1)
          lt     = indmtr(ideb,j,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_2d.for"
#include "FastS/Compute/LU/lu_dinv_2d.for"
        enddo

      !!! Fin plan kdeb

      do j= jdeb+ipas,jfin,ipas
      do i= ideb+ipas,ifin,ipas

          l      = inddm(i,j,1)
          lt     = indmtr(i,j,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_2d.for"
#include "FastS/Compute/LU/lu_j_2d.for"
#include "FastS/Compute/LU/lu_dinv_2d.for"
        enddo
        enddo

      ENDIF
 
      end
