c***********************************************************************
c     $Date: 2010-07-12 18:57:35 +0200 (Mon, 12 Jul 2010) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlu_ale_u_SA(ndom,  param_int, param_real,
     &                      visco,sa_real,
     &                      ind_loop,
     &                      drodm_out,rop,rop_1,
     &                      ti,tj,tk,venti,ventj,ventk,
     &                      coe, mjrnewton)
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

#include "FastS/param_solver.h"

      INTEGER_E ndom, ind_loop(6), param_int(0:*), mjrnewton
 

      REAL_E visco(5),sa_real(3)

      REAL_E  param_real(0:*)
      REAL_E drodm_out(param_int(NDIMDX)  ,param_int(NEQ)),
     &       coe(param_int(NDIMDX)  ,param_int(NEQ_COE)),
     &       rop(param_int(NDIMDX)  ,param_int(NEQ)),
     &       rop_1(param_int(NDIMDX),param_int(NEQ))
      REAL_E ti(param_int(NDIMDX_MTR),param_int(NEQ_IJ)),
     &       tj(param_int(NDIMDX_MTR),param_int(NEQ_IJ)),
     &       tk(param_int(NDIMDX_MTR),param_int(NEQ_K))

      REAL_E venti(param_int(NDIMDX_VENT)*param_int(NEQ_VENT))
      REAL_E ventj(param_int(NDIMDX_VENT)*param_int(NEQ_VENT))
      REAL_E ventk(param_int(NDIMDX_VENT)*param_int(NEQ_VENT))

c Var loc
      INTEGER_E  inci,incj,inck,l,i,j,k,kdmax,kd,lmax,ll,ndo,
     & kddeb,kdfin,ipas,kfin,kdeb,jfin,jdeb,ifin,ideb,
     & l1,l2,lt,lt1,lt2,lij,lv,lv1,lv2,v2ven,v3ven,
     & inci_ven,incj_ven,inck_ven,inci2_ven,incj2_ven,inck2_ven,
     & inci2_mtr,incj2_mtr,inck2_mtr,inci_mtr,incj_mtr,inck_mtr,
     & ls,l1s,incis,incjs,incks

      REAL_E gam2,gam1,gamm1,cp,xal,diag, ratiom,
     & b11,b12,b13,b14,b15,b21,b22,b23,b24,b25,b31,b32,b33,b34,b35,b41,
     & b42,b43,b44,b45,b51,b52,b53,b54,b55,
     & b1,b2,b3,b4,b5,b6,ck_vent,ue,ve,we,qen,
     & signe,r,u,v,w,t,q2,h,ph2,qn,tcx,tcy,tcz,
     & anulam,temp01,cmus1,coesut,
     & ro_old,u_old,v_old,w_old,t_old,roe_old,nu_old,r_1,cvinv,cvinv2

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

      gam1    = param_real(GAMMA)
      gamm1   = gam1 - 1.
      cp      = param_real(GAMMA)*param_real(CVINF)
      gam2    = gamm1-1.

      cvinv   = 1./param_real(CVINF)
      cvinv2  = 0.5*cvinv

      cmus1  = visco(5)
      temp01 = 1./visco(4)
      coesut = visco(3) * (1.+cmus1*temp01)
      ratiom = sa_real(SA_RATIOM)
 
       inci      =-1
       incj      =-param_int(NIJK)
       inck      =-param_int(NIJK)*param_int(NIJK+1)
       inci_mtr  =-param_int(NIJK_MTR)
       incj_mtr  =-param_int(NIJK_MTR+1)
       inck_mtr  =-param_int(NIJK_MTR+2)

       inci_ven  =-param_int(NIJK_VENT)
       incj_ven  =-param_int(NIJK_VENT+1)
       inck_ven  =-param_int(NIJK_VENT+2)

       inci2_mtr = 2*inci_mtr
       incj2_mtr = 2*incj_mtr
       inck2_mtr = 2*inck_mtr
       inci2_ven = 2*inci_ven
       incj2_ven = 2*incj_ven
       inck2_ven = 2*inck_ven

       ipas      =-1
       signe     =-0.5
      
       if(param_int(NEQ_VENT).eq.2) then
         ck_vent =0.
       else
         ck_vent =1.
       endif

       v2ven =   param_int(NDIMDX_VENT)
       v3ven = 2*param_int(NDIMDX_VENT)*ck_vent

        kfin  = ind_loop(5)
        jfin  = ind_loop(3)
        ifin  = ind_loop(1)
        kdeb  = ind_loop(6)
        jdeb  = ind_loop(4)
        ideb  = ind_loop(2)

        incis = inci
        incjs = incj
        incks = inck

      IF(param_int(ITYPZONE).eq.0) THEN !domaine 3d general

       !!Diag
       do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, kdeb)
             ls = l
#include     "FastS/Compute/LU/lu_d_SA.for"
           enddo
        enddo

      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)
        l = inddm(ideb,jdeb,kdeb)
        ls = l
#include "FastS/Compute/LU/lu_dinv_SA.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          ls = l
          lt     = indmtr(i,jdeb,kdeb)
          lv     = indven(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include "FastS/Compute/LU/lu_dinv_SA.for"
        enddo

        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =  inddm(ideb,j,kdeb)
          ls = l
          lt     = indmtr(ideb,j,kdeb)
          lv     = indven(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include "FastS/Compute/LU/lu_dinv_SA.for"

          do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,kdeb)
            ls = l
            lt     = indmtr(i,j,kdeb)
            lv     = indven(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_j_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

          enddo
        enddo
        !!! plan kdeb termine

        if (mjrnewton == 1) then
!!mise a jour Newton 1er plan k
           do j= jdeb,jfin,ipas
C     DIR$ IVDEP
              do  i= ideb,ifin,ipas

                 l  = inddm(i, j, kdeb)
                 ls = l
#include     "FastS/Compute/LU/mjr_newton_SA.for"
              enddo
           enddo
        endif

        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, k)
             ls = l
#include     "FastS/Compute/LU/lu_d_SA.for"
           enddo
          enddo

          l      =  inddm(ideb,jdeb,k)
          ls = l
          lt     = indmtr(ideb,jdeb,k)
          lv     = indven(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =  inddm(i,jdeb,k)
             ls = l
             lt     = indmtr(i,jdeb,k)
             lv     = indven(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_k_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =  inddm(ideb,j,k)
             ls = l
             lt     = indmtr(ideb,j,k)
             lv     = indven(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_k_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

             do  i= ideb+ipas,ifin,ipas
      
               l = inddm(i,j,k)
               ls = l
               lt= indmtr(i,j,k)
               lv= indven(i,j,k)

                xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_j_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_k_3dfull_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_dinv_SA.for"
             enddo
          enddo

          if (mjrnewton == 1) then
             do j= jfin,jdeb
C     DIR$ IVDEP
                do  i= ifin, ideb

                   l  = inddm(i, j, k)
                   ls = l
#include     "FastS/Compute/LU/mjr_newton_SA.for"
                enddo
             enddo
          endif

        enddo


      ELSEIF(param_int(ITYPZONE).eq.1) THEN !maillage 3d k homogene_SA:

       !!Diag
       do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, kdeb)
             ls = l
#include     "FastS/Compute/LU/lu_d_SA.for"
           enddo
        enddo

      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)
        l = inddm(ideb,jdeb,kdeb)
        ls = l
#include "FastS/Compute/LU/lu_dinv_SA.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          ls = l
          lt     = indmtr(i,jdeb,kdeb)
          lv     = indven(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include "FastS/Compute/LU/lu_dinv_SA.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =  inddm(ideb,j,kdeb)
          ls = l
          lt     = indmtr(ideb,j,kdeb)
          lv     = indven(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include "FastS/Compute/LU/lu_dinv_SA.for"

          do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,kdeb)
            ls = l
            lt     = indmtr(i,j,kdeb)
            lv     = indven(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_j_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

          enddo
        enddo
        !!! plan kdeb termine

        if (mjrnewton == 1) then
!!mise a jour Newton 1er plan k
           do j= jdeb,jfin,ipas
C     DIR$ IVDEP
              do  i= ideb,ifin,ipas

                 l  = inddm(i, j, kdeb)
                 ls = l
#include     "FastS/Compute/LU/mjr_newton_SA.for"
              enddo
           enddo
        endif

        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, k)
             ls = l
#include     "FastS/Compute/LU/lu_d_SA.for"
           enddo
          enddo

          l      =  inddm(ideb,jdeb,k)
          ls = l
          lt     = indmtr(ideb,jdeb,k)
          lv     = indven(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =  inddm(i,jdeb,k)
             ls = l
             lt     = indmtr(i,jdeb,k)
             lv     = indven(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_k_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =  inddm(ideb,j,k)
             ls = l
             lt     = indmtr(ideb,j,k)
             lv     = indven(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_k_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

             do  i= ideb+ipas,ifin,ipas
      
               l = inddm(i,j,k)
               ls = l
               lt= indmtr(i,j,k)
               lv= indven(i,j,k)

                xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_j_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_k_3dhomogene_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_dinv_SA.for"
             enddo
          enddo

          if (mjrnewton == 1) then
             do j= jfin,jdeb
C     DIR$ IVDEP
                do  i= ifin, ideb

                   l  = inddm(i, j, k)
                   ls = l
#include     "FastS/Compute/LU/mjr_newton_SA.for"
                enddo
             enddo
          endif

        enddo

      ELSEIF(param_int(ITYPZONE).eq.2) THEN !maillage 3d cart_SAesien



       tcx = ti(1,1)
       tcy = tj(1,1)
       tcz = tk(1,1)

       !!Diag
       do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, kdeb)
             ls = l
#include     "FastS/Compute/LU/lu_d_SA.for"
           enddo
        enddo

      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)
        l = inddm(ideb,jdeb,kdeb)
        ls = l
#include "FastS/Compute/LU/lu_dinv_SA.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,kdeb)
          ls = l
          lt     = indmtr(i,jdeb,kdeb)
          lv     = indven(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include "FastS/Compute/LU/lu_dinv_SA.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =  inddm(ideb,j,kdeb)
          ls = l
          lt     = indmtr(ideb,j,kdeb)
          lv     = indven(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include "FastS/Compute/LU/lu_dinv_SA.for"

          do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,kdeb)
            ls = l
            lt     = indmtr(i,j,kdeb)
            lv     = indven(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_j_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

          enddo
        enddo
        !!! plan kdeb termine

        if (mjrnewton == 1) then
!!mise a jour Newton 1er plan k
           do j= jdeb,jfin,ipas
C     DIR$ IVDEP
              do  i= ideb,ifin,ipas

                 l  = inddm(i, j, kdeb)
                 ls = l
#include     "FastS/Compute/LU/mjr_newton_SA.for"
              enddo
           enddo
        endif

        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, k)
             ls = l
#include     "FastS/Compute/LU/lu_d_SA.for"
           enddo
          enddo

          l      =  inddm(ideb,jdeb,k)
          ls = l
          lt     = indmtr(ideb,jdeb,k)
          lv     = indven(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =  inddm(i,jdeb,k)
             ls = l
             lt     = indmtr(i,jdeb,k)
             lv     = indven(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_k_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =  inddm(ideb,j,k)
             ls = l
             lt     = indmtr(ideb,j,k)
             lv     = indven(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_k_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include    "FastS/Compute/LU/lu_dinv_SA.for"

             do  i= ideb+ipas,ifin,ipas
      
               l = inddm(i,j,k)
               ls = l
               lt= indmtr(i,j,k)
               lv= indven(i,j,k)

                xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_j_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_k_3dcart_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_SA.for"
#include       "FastS/Compute/LU/lu_dinv_SA.for"
             enddo
          enddo
   
          if (mjrnewton == 1) then
!! mise a jour newton
             do j= jfin,jdeb
C     DIR$ IVDEP
                do  i= ifin, ideb

                   l  = inddm(i, j, k)
                   ls = l
#include     "FastS/Compute/LU/mjr_newton_SA.for"
                enddo
             enddo
          endif

        enddo


      ELSE !2d


       !!Diag
       do j= jdeb,jfin,ipas
           do  i= ideb,ifin,ipas

             l  = inddm(i, j, kdeb)
             ls     = l
#include     "FastS/Compute/LU/lu_d_2d_SA.for"
           enddo
       enddo

       !!! coin (ideb,jdeb,kdb)
       l = inddm(ideb,jdeb,1)
       ls     = l
#include "FastS/Compute/LU/lu_dinv_2d_SA.for"

      !!! ligne (jdeb,kdeb) dans le plan kdeb
       do i= ideb+ipas,ifin,ipas

          l      = inddm(i,jdeb,1)
          ls     = l
          lt     = indmtr(i,jdeb,1)
          lv     = indven(i,jdeb,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_2d_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_2d_SA.for"
#include "FastS/Compute/LU/lu_dinv_2d_SA.for"
       enddo

       do j= jdeb+ipas,jfin,ipas

          l      =  inddm(ideb,j,1)
          ls     = l
          lt     = indmtr(ideb,j,1)
          lv     = indven(ideb,j,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_2d_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_2d_SA.for"
#include "FastS/Compute/LU/lu_dinv_2d_SA.for"
       enddo
       do j= jdeb+ipas,jfin,ipas
         do i= ideb+ipas,ifin,ipas

            l      =  inddm(i,j,1)
            ls     = l
            lt     = indmtr(i,j,1)
            lv     = indven(i,j,1)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_2d_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_2d_SA.for"
#include    "FastS/Compute/LU/lu_j_2d_SA_ale.for"
#include "FastS/Compute/LU/mjr_drodm_2d_SA.for"
#include    "FastS/Compute/LU/lu_dinv_2d_SA.for"
          enddo
       enddo

       if (mjrnewton == 1) then
!! mise a jour newton
          do j= jdeb,jfin,ipas
C     DIR$ IVDEP
             do  i= ifin, ideb

                l  = inddm(i, j, 1)
                ls     = l
#include     "FastS/Compute/LU/mjr_newton_2d_SA.for"
             enddo
          enddo
       endif


      ENDIF
 
      end
