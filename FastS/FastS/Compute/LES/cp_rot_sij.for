c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 févr. 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_rot_sij(ndom, param_int, neq_rot,
     &                     ind_loop, 
     &                     xmut,rop, ti, tj, tk, vol, rot)
c***********************************************************************
c_P                          O N E R A
c     ACT
c_A    Calcul de la contribution des termes sources 
c_A    pour l'equation de Spalart Allmaras 
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    Navier-Stokes
c
c     INP
c_I    ndom      : numero du domaine calcule
c_I    dvardc(6) : gradients de nutild primitif aux centres
c_I                des cellules
c_I    dvardc(5) : gradients de ronutild conservatif aux centres
c_I                des cellules
c_I    vol       : volumes
c_I    rotn      : norme du rotationnel
c_I    dlng      : distance a la paroi
c
c     OUT
c
c     I/O
c_I    drodm    : terme source de l'equation de Spalart Allmaras 
c***********************************************************************
      implicit  none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), neq_rot, param_int(0:*)

      REAL_E xmut( param_int(NDIMDX) )
      REAL_E rot ( param_int(NDIMDX) * neq_rot)
      REAL_E rop( param_int(NDIMDX)  * param_int(NEQ) )

      REAL_E ti( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) * param_int(NEQ_K ) )

      REAL_E vol(param_int(NDIMDX_MTR))

c Var loc 
      INTEGER_E incmax,i,j,k,l,icoe_pos,inci,incj,inck,inci_mtr,
     & incj_mtr,inck_mtr, l1,l2,l3,l4,l5,l6,ltij,lij,lt,lvo,
     & v1,v2,v3,v4,v1mtr,v2mtr,v3mtr

      REAL_E  tjx,tjy,tjz,tjx1,tjy1,tjz1,
     & tix,tiy,tiz,tix1,tiy1,tiz1,
     & tkx,tky,tkz,tkx1,tky1,tkz1,u1,u2,u3,u4,u5,u6,xvol,
     & dudx,dudy,dudz,
     & dvdx,dvdy,dvdz,
     & dwdx,dwdy,dwdz


#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      incmax=param_int(NIJK)*param_int(NIJK+1)*param_int(NIJK+4)

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)

      v1 =   0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)

      v1mtr =   0
      v2mtr =   param_int(NDIMDX_MTR)
      v3mtr = 2*param_int(NDIMDX_MTR)

      IF(param_int(ITYPZONE).eq.0) THEN  !domaine 3d:

#include   "FastS/Compute/loop_begin.for"
#include       "FastS/Compute/LES/metric_3dfull.for"
#include       "FastS/Compute/LES/rot_3dfull.for" 
#include   "FastS/Compute/loop_end.for"

      !!!!!
      !!!!!
      !!!!!
      !3dhomo
      ELSEIF(param_int(ITYPZONE).eq.1)  THEN 

#include   "FastS/Compute/loop_begin.for"
#include       "FastS/Compute/LES/metric_3dhomo.for"
#include       "FastS/Compute/LES/rot_3dhomo.for" 
#include   "FastS/Compute/loop_end.for"

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
          tix = ti(lt)
          tjy = tj(lt)
          tkz = tk(lt)
          xvol    = 0.5/vol(lvo)


           do k = ind_loop(5), ind_loop(6)
            do j = ind_loop(3), ind_loop(4)
             lij  =       inddm( ind_loop(1) , j, k)
!DEC$ IVDEP
             do l = lij, lij +  ind_loop(2)- ind_loop(1)
#include       "FastS/Compute/LES/rot_3dcart.for" 
             enddo
            enddo
           enddo
      !!!!!
      !!!!!
      !!!!!
      !!!!!
      !2D     
      ELSE !2d:  pas de 2D en LES

!$OMP SINGLE
           write(*,*)'calcul LES et zone 2D...Probleme'
!$OMP END SINGLE

      ENDIF! typezone (2d,3d, cart,...0

      end
