c***********************************************************************
c     $Date: 2010-07-12 18:57:35 +0200 (Mon, 12 Jul 2010) $
c     $Revision: 40 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine invlussor_ale_l(ndom, param_int, param_real,
     &                   ind_loop, ind_loop_sdm,
     &                   drodm_out,rop,
     &                   ti,tj,tk, venti,ventj,ventk,
     &                   coe,ssor, ssor_size)
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

      INTEGER_E ndom, ind_loop(6), param_int(0:*), ssor_size,
     &     ind_loop_sdm(6)
 
      REAL_E  param_real(0:*)
      REAL_E drodm_out(ssor_size,param_int(NEQ)),
     &       coe(param_int(NDIMDX),param_int(NEQ_COE)),
     &       rop(param_int(NDIMDX),param_int(NEQ)),
     &     ssor(ssor_size,param_int(NEQ))

      REAL_E ti(param_int(NDIMDX_MTR),param_int(NEQ_IJ)),
     &       tj(param_int(NDIMDX_MTR),param_int(NEQ_IJ)),
     &       tk(param_int(NDIMDX_MTR),param_int(NEQ_K))

      REAL_E venti(param_int(NDIMDX_VENT)*param_int(NEQ_VENT))
      REAL_E ventj(param_int(NDIMDX_VENT)*param_int(NEQ_VENT))
      REAL_E ventk(param_int(NDIMDX_VENT)*param_int(NEQ_VENT))

c Var loc
      INTEGER_E  inci,incj,inck,l,i,j,k,kdmax,kd,lmax,ll,ndo,
     & kddeb,kdfin,ipas,kfin,kdeb,jfin,jdeb,ifin,ideb,
     & l1,l2,lt,lt1,lt2,lv,lv1,lv2,v2ven,v3ven,
     & inci_ven,incj_ven,inck_ven,inci2_ven,incj2_ven,inck2_ven,
     & inci2_mtr,incj2_mtr,inck2_mtr,inci_mtr,incj_mtr,inck_mtr,
     & ls,l1s,incis,incjs,incks

      REAL_E gam2,gam1,gamm1,cp,xal,diag,
     & b11,b12,b13,b14,b15,b21,b22,b23,b24,b25,b31,b32,b33,b34,b35,b41,
     & b42,b43,b44,b45,b51,b52,b53,b54,b55,
     & b1,b2,b3,b4,b5,ck_vent,ue,ve,we,qen,
     & signe,r,u,v,w,t,q2,h,ph2,qn,tcx,tcy,tcz

#include "FastS/formule_param.h"
#include "FastS/formule_ssor_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

      gam1    = param_real(GAMMA)
      gamm1   = gam1 - 1.
      cp      = param_real(GAMMA)*param_real(CVINF)
      gam2    = gamm1-1.
 
      inci      = 1
      incj      = param_int(NIJK)
      inck      = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr  = param_int(NIJK_MTR)
      incj_mtr  = param_int(NIJK_MTR+1)
      inck_mtr  = param_int(NIJK_MTR+2)
      inci2_mtr = 0
      incj2_mtr = 0
      inck2_mtr = 0

      inci_ven  = param_int(NIJK_VENT)
      incj_ven  = param_int(NIJK_VENT+1)
      inck_ven  = param_int(NIJK_VENT+2)
      inci2_ven = 0
      incj2_ven = 0
      inck2_ven = 0

      ipas      = 1
      signe     = 0.5
      
      if(param_int(NEQ_VENT).eq.2) then
         ck_vent =0.
      else
         ck_vent =1.
      endif

      v2ven =   param_int(NDIMDX_VENT)
      v3ven = 2*param_int(NDIMDX_VENT)*ck_vent


      kdeb  = ind_loop(5)
      jdeb  = ind_loop(3)
      ideb  = ind_loop(1)
      kfin  = ind_loop(6)
      jfin  = ind_loop(4)
      ifin  = ind_loop(2)

      ssor_i= ind_loop_sdm(2)-ind_loop_sdm(1)+1 +2*param_int(NIJK+3) !taille de la fenetre + ghostcells
      ssor_j= ind_loop_sdm(4)-ind_loop_sdm(3)+1 +2*param_int(NIJK+3)

      incis = 1
      incjs = ssor_i
      incks = ssor_i * ssor_j

      IF(param_int(ITYPZONE).eq.0) THEN !domaine 3d general


      !!! coin (ideb,jdeb,kdb)
        l =   inddm(ideb,jdeb,kdeb)
        ls= indssor(ideb,jdeb,kdeb)

#include "FastS/Compute/LU/lu_dinv.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      =   inddm(i,jdeb,kdeb)
          lt     =  indmtr(i,jdeb,kdeb)
          lv     =  indven(i,jdeb,kdeb)
          ls     = indssor(i,jdeb,kdeb)
  
          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dfull_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =   inddm(ideb,j,kdeb)
          lt     =  indmtr(ideb,j,kdeb)
          lv     =  indven(ideb,j,kdeb)
          ls     = indssor(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dfull_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"

          do i= ideb+ipas,ifin,ipas

            l      =   inddm(i,j,kdeb)
            lt     =  indmtr(i,j,kdeb)
            lv     =  indven(i,j,kdeb)
            ls     = indssor(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dfull_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_j_3dfull_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
        enddo
        !!! plan kdeb termine
        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          l      =   inddm(ideb,jdeb,k)
          lt     =  indmtr(ideb,jdeb,k)
          lv     =  indven(ideb,jdeb,k)
          ls     = indssor(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dfull_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =   inddm(i,jdeb,k)
             lt     =  indmtr(i,jdeb,k)
             lv     =  indven(i,jdeb,k)
             ls     = indssor(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dfull_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_k_3dfull_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =   inddm(ideb,j,k)
             lt     =  indmtr(ideb,j,k)
             lv     =  indven(ideb,j,k)
             ls     = indssor(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dfull_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_k_3dfull_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"

             do  i= ideb+ipas,ifin,ipas
      
               l      =   inddm(i,j,k)
               lt     =  indmtr(i,j,k)
               lv     =  indven(i,j,k)
               ls     = indssor(i,j,k)

               xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dfull_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_j_3dfull_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_k_3dfull_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_dinv.for"
             enddo
          enddo
        enddo


      ELSEIF(param_int(ITYPZONE).eq.1) THEN !maillage 3d k homogene:

      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)
        l =   inddm(ideb,jdeb,kdeb)
        ls= indssor(ideb,jdeb,kdeb)

#include "FastS/Compute/LU/lu_dinv.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      =   inddm(i,jdeb,kdeb)
          lt     =  indmtr(i,jdeb,kdeb)
          lv     =  indven(i,jdeb,kdeb)
          ls     = indssor(i,jdeb,kdeb)
  
          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dhomogene_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =   inddm(ideb,j,kdeb)
          lt     =  indmtr(ideb,j,kdeb)
          lv     =  indven(ideb,j,kdeb)
          ls     = indssor(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dhomogene_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"

          do i= ideb+ipas,ifin,ipas

            l      =   inddm(i,j,kdeb)
            lt     =  indmtr(i,j,kdeb)
            lv     =  indven(i,j,kdeb)
            ls     = indssor(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dhomogene_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_j_3dhomogene_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
        enddo
        !!! plan kdeb termine
        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          l      =   inddm(ideb,jdeb,k)
          lt     =  indmtr(ideb,jdeb,k)
          lv     =  indven(ideb,jdeb,k)
          ls     = indssor(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dhomogene_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =   inddm(i,jdeb,k)
             lt     =  indmtr(i,jdeb,k)
             lv     =  indven(i,jdeb,k)
             ls     = indssor(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dhomogene_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_k_3dhomogene_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =   inddm(ideb,j,k)
             lt     =  indmtr(ideb,j,k)
             lv     =  indven(ideb,j,k)
             ls     = indssor(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dhomogene_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_k_3dhomogene_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"

             do  i= ideb+ipas,ifin,ipas
      
               l      =   inddm(i,j,k)
               lt     =  indmtr(i,j,k)
               lv     =  indven(i,j,k)
               ls     = indssor(i,j,k)

                xal   = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dhomogene_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_j_3dhomogene_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_k_3dhomogene_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_dinv.for"
             enddo
          enddo
        enddo

      ELSEIF(param_int(ITYPZONE).eq.2) THEN !maillage 3d cartesien


       tcx = ti(1,1)
       tcy = tj(1,1)
       tcz = tk(1,1)
      !!! on parcourt le domaine en 7 passes pour traiter les bord sans mettre a zero le drodm sur maille fictive

      !!! coin (ideb,jdeb,kdb)
        l =   inddm(ideb,jdeb,kdeb)
        ls= indssor(ideb,jdeb,kdeb)

#include "FastS/Compute/LU/lu_dinv.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      =   inddm(i,jdeb,kdeb)
          lv     =  indven(i,jdeb,kdeb)
          ls     = indssor(i,jdeb,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_3dcart_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"
        enddo


        do j= jdeb+ipas,jfin,ipas


          !!! ligne (ideb,kdeb) dans le plan kdeb
          l      =   inddm(ideb,j,kdeb)
          lv     =  indven(ideb,j,kdeb)
          ls     = indssor(ideb,j,kdeb)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_3dcart_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"

          do i= ideb+ipas,ifin,ipas

            l      =   inddm(i,j,kdeb)
            lv     =  indven(i,j,kdeb)
            ls     = indssor(i,j,kdeb)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dcart_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_j_3dcart_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
        enddo
        !!! plan kdeb termine
        !!! le domaine sans les mailles du bord
        do  k= kdeb+ipas,kfin,ipas

          l      =   inddm(ideb,jdeb,k)
          lv     =  indven(ideb,jdeb,k)
          ls     = indssor(ideb,jdeb,k)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_k_3dcart_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor.for"
#include "FastS/Compute/LU/lu_dinv.for"

          !!! Fin plan jdeb 
          do i= ideb+ipas,ifin,ipas

             l      =   inddm(i,jdeb,k)
             lv     =  indven(i,jdeb,k)
             ls     = indssor(i,jdeb,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_3dcart_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_k_3dcart_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"
          enddo
          do  j= jdeb+ipas,jfin,ipas

             l      =   inddm(ideb,j,k)
             lv     =  indven(ideb,j,k)
             ls     = indssor(ideb,j,k)

             xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_j_3dcart_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_k_3dcart_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor.for"
#include    "FastS/Compute/LU/lu_dinv.for"

             do  i= ideb+ipas,ifin,ipas
      
               l      =   inddm(i,j,k)
               lv     =  indven(i,j,k)
               ls     = indssor(i,j,k)

               xal    = coe(l,1)*signe

#include       "FastS/Compute/LU/lu_i_3dcart_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_j_3dcart_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_k_3dcart_ale.for"
#include       "FastS/Compute/LU/mjr_drodm_ssor.for"
#include       "FastS/Compute/LU/lu_dinv.for"
             enddo
          enddo
        enddo



      ELSE !2d

        !!! coin (ideb,jdeb,kdb)
        l =   inddm(ideb,jdeb,1)
        ls= indssor(ideb,jdeb,1)
#include "FastS/Compute/LU/lu_dinv_2d.for"

        !!! ligne (jdeb,kdeb) dans le plan kdeb
        do i= ideb+ipas,ifin,ipas

          l      =   inddm(i,jdeb,1)
          lt     =  indmtr(i,jdeb,1)
          lv     =  indven(i,jdeb,1)
          ls     = indssor(i,jdeb,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_i_2d_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor_2d.for"
#include "FastS/Compute/LU/lu_dinv_2d.for"
        enddo

        do j= jdeb+ipas,jfin,ipas

          l      =   inddm(ideb,j,1)
          lt     =  indmtr(ideb,j,1)
          lv     =  indven(ideb,j,1)
          ls     = indssor(ideb,j,1)

          xal    = coe(l,1)*signe

#include "FastS/Compute/LU/lu_j_2d_ale.for"
#include "FastS/Compute/LU/mjr_drodm_ssor_2d.for"
#include "FastS/Compute/LU/lu_dinv_2d.for"
          do i= ideb+ipas,ifin,ipas

            l      =   inddm(i,j,1)
            lt     =  indmtr(i,j,1)
            lv     =  indven(i,j,1)
            ls     = indssor(i,j,1)

            xal    = coe(l,1)*signe

#include    "FastS/Compute/LU/lu_i_2d_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor_2d.for"
#include    "FastS/Compute/LU/lu_j_2d_ale.for"
#include    "FastS/Compute/LU/mjr_drodm_ssor_2d.for"
#include    "FastS/Compute/LU/lu_dinv_2d.for"
          enddo
        enddo

      ENDIF
 
      end
