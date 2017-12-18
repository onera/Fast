c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 févr. 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine spsource_MEANFLOW(ndom, nitcfg, param_int, param_real,
     &                     ind_loop, rop, coe, vol, drodm)
c***********************************************************************
c_P                          O N E R A
c     ACT
c_A    Calcul de la contribution du terme source centrifuge 
c_A    pour les equations du champ moyen (moment)
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    Navier-Stokes
c
c     INP
c_I    ndom      : numero du domaine calcule
c_I    vol       : volumes

c
c     OUT
c
c     I/O
c_I    drodm    : terme source champ moyen (moment)
c***********************************************************************
      implicit  none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), nitcfg, param_int(0:*)

      REAL_E rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E coe( param_int(NDIMDX) , param_int(NEQ_COE) )
      REAL_E drodm( param_int(NDIMDX), param_int(NEQ) )
      REAL_E vol(param_int(NDIMDX_MTR))

      REAL_E param_real(0:*)
      
c Var loc 
      INTEGER_E incmax,i,j,k,l,icoe_pos,inci,incj,inck,inci_mtr,
     & incj_mtr,inck_mtr, l1,l2,l3,l4,l5,l6,ltij,lij,lt,ndimdx,lvo

      REAL_E axe(3),tsource(3)
      REAL_E omg,omg1,omg2,omg3,coe2,coe3,coe4
      
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      incmax=param_int(NIJK)*param_int(NIJK+1)*param_int(NIJK+4)

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)

      ndimdx   = param_int(NDIMDX)

      ! Donnees cinematiques du terme source 
      axe(1) = param_real(ROT_OMEGA)
      axe(2) = param_real(ROT_OMEGA+1)
      axe(3) = param_real(ROT_OMEGA+2)
      
      omg = param_real(ROT_TETAP)
      
      IF (axe(2).eq.0. .and. axe(3).eq.0.) THEN
        omg1 = omg
        omg2 = 0.
        omg3 = 0.
      ELSE IF (axe(1).eq.0. .and. axe(3).eq.0.) THEN
        omg1 = 0.
        omg2 = omg
        omg3 = 0.
      ELSE IF (axe(1).eq.0. .and. axe(2).eq.0.) THEN
        omg1 = 0.
        omg2 = 0.
        omg3 = omg
      ENDIF
     
      IF(param_int(ITYPZONE).eq.0) THEN !domaine 3d:

         If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor

#include  "FastS/Compute/loop_begin.for"

#include       "FastS/Compute/ALE/source_centrifuge.for"
#include       "FastS/Compute/ALE/source_centrifuge_LU.for"
            drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
            drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
            drodm(l,4)= drodm(l,4) + vol(lvo)*tsource(3)
            
#include  "FastS/Compute/loop_end.for"

         Else                   !calcul explicit, Stockage terme source inutile
            
#include  "FastS/Compute/loop_begin.for"
            
#include       "FastS/Compute/ALE/source_centrifuge.for"
            drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
            drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
            drodm(l,4)= drodm(l,4) + vol(lvo)*tsource(3)
            
#include  "FastS/Compute/loop_end.for"
            
            
         Endif                  !explicite/implicite
       
!!!!!
!!!!!
!!!!!
!!!!!
!!!!!
!3dhomo
      ELSEIF(param_int(ITYPZONE).eq.1)  THEN 
         
         
         If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor

#include  "FastS/Compute/loop_begin.for"
            
#include       "FastS/Compute/ALE/source_centrifuge.for"
#include       "FastS/Compute/ALE/source_centrifuge_LU.for"
            drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
            drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
            drodm(l,4)= drodm(l,4) + vol(lvo)*tsource(3)
#include  "FastS/Compute/loop_end.for"
            
         Else                   !calcul explicit, Stockage terme source inutile
            
#include  "FastS/Compute/loop_begin.for"
            
#include       "FastS/Compute/ALE/source_centrifuge.for"
            drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
            drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
            drodm(l,4)= drodm(l,4) + vol(lvo)*tsource(3)
#include  "FastS/Compute/loop_end.for"
            

         Endif                  !explicite/implicite 3dhomo
         
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !3dcart
      ELSEIF(param_int(ITYPZONE).eq.2)  THEN 
         

      !metric
         lt  = indmtr(1 , 1, 1)
         lvo = lt
         
         If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor
            
            do k = ind_loop(5), ind_loop(6)
               do j = ind_loop(3), ind_loop(4)
                  lij  =       inddm( ind_loop(1) , j, k)
!     DEC$ IVDEP
                  do l = lij, lij +  ind_loop(2)- ind_loop(1)
                     
#include       "FastS/Compute/ALE/source_centrifuge.for"
#include       "FastS/Compute/ALE/source_centrifuge_LU.for"
                     drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
                     drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
                     drodm(l,4)= drodm(l,4) + vol(lvo)*tsource(3)
                  enddo
               enddo
            enddo

         Else                   !calcul explicit, Stockage terme source inutile
            
            do k = ind_loop(5), ind_loop(6)
               do j = ind_loop(3), ind_loop(4)
                  lij  =       inddm( ind_loop(1) , j, k)
!     DEC$ IVDEP
                  do l = lij, lij +  ind_loop(2)- ind_loop(1)
                     
#include       "FastS/Compute/ALE/source_centrifuge.for"
                     drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
                     drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
                     drodm(l,4)= drodm(l,4) + vol(lvo)*tsource(3)
                  enddo
               enddo
            enddo

         Endif                  !explicite/implicite 3dcart
         
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !2D     
      ELSE !2d:

       If(param_int(ITYPCP).le.1.and.nitcfg.le.1) then !calcul implicite, on stocke coe pour ssor

          do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
             do l = lij, lij +  ind_loop(2)- ind_loop(1)

#include       "FastS/Compute/ALE/source_centrifuge_2d.for"
#include       "FastS/Compute/ALE/source_centrifuge_LU.for"
               drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
               drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
           enddo
           enddo
          enddo

       Else  !calcul explicit, Stockage terme source inutile


          do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
             do l = lij, lij +  ind_loop(2)- ind_loop(1)

#include       "FastS/Compute/ALE/source_centrifuge_2d.for"
               drodm(l,2)= drodm(l,2) + vol(lvo)*tsource(1)
               drodm(l,3)= drodm(l,3) + vol(lvo)*tsource(2)
           enddo
           enddo
          enddo
 
       endif!explicite/implicite 2d
      ENDIF! typezone (2d,3d, cart,...0
      end
