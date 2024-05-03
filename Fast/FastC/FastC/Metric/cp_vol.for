c    ***********************************************************************
c     $Date: 2010-06-14 09:57:46 +0200 (Mon, 14 Jun 2010) $
c     $Revision: 60 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_vol( param_int,
     &                    x, y, z, ti, tj, tk, ti0, tj0, tk0, vol,
     &                    ind_loop )
c***********************************************************************
c_P                          O N E R A
c
c     ACT
c_A    calcul metrique: normale + surface + volume en volume finis
c
c     VAL
c
c     INP
C_I    neq_ij,neq_k : dimension metrique
c_/    x,y,z        : coordonnees du maillage
c
c     OUT
c     ti,tj,tk,si,sj,sk,vol
c***********************************************************************
      implicit none

#include "FastC/param_solver.h"

      INTEGER_E ind_loop(6),param_int(0:*)

      REAL_E x(param_int( NDIMDX_XYZ )),y(param_int( NDIMDX_XYZ )),
     &       z(param_int( NDIMDX_XYZ ))    
C_OUT
      REAL_E  ti(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E ti0(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E  tj(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E tj0(param_int( NDIMDX_MTR ),param_int( NEQ_IJ )) 
      REAL_E  tk(param_int( NDIMDX_MTR ),param_int( NEQ_K  )) 
      REAL_E tk0(param_int( NDIMDX_MTR ),param_int( NEQ_K  )) 
      REAL_E vol(param_int( NDIMDX_MTR )) 

c Var loc 
      INTEGER_E k,j,i,l0,l1,l5,l4,l2,l3,l7,l6,incmax,
     & lip,ljp,lkp,l,inci,incj,inck,ne,i1,i2,j1,j2,k1,k2,lmax,
     & inci_mtr,incj_mtr,inck_mtr,iv,jv,kv

      REAL_E ax,ay,az,tk1,tk2,dz,
     & tix10,tix11,tiy10,tiy11,tiz10,tiz11,
     & tix20,tix21,tiy20,tiy21,tiz20,tiz21,
     & tjx10,tjx11,tjy10,tjy11,tjz10,tjz11,
     & tjx20,tjx21,tjy20,tjy21,tjz20,tjz21,
     & tkx10,tkx11,tky10,tky11,tkz10,tkz11,
     & tkx20,tkx21,tky20,tky21,tkz20,tkz21


#include "FastC/formule_mtr_param.h"
#include "FastC/formule_xyz_param.h"

      inci = 1
      incj = param_int( NIJK_XYZ ) 
      inck = param_int( NIJK_XYZ ) *param_int( NIJK_XYZ +1 ) 

      iv = param_int( IJKV   )
      jv = param_int( IJKV +1)
      kv = param_int( IJKV +2)

      inci_mtr = param_int( NIJK_MTR   )
      incj_mtr = param_int( NIJK_MTR +1)
      inck_mtr = param_int( NIJK_MTR +2)

c-----le nbr de metrique varie selon le type de domaine
      IF(param_int( ITYPZONE ).eq.0) THEN !Domaine 3D quelconque

        !calcul volume
        do k=ind_loop(5),ind_loop(6)
        do j=ind_loop(3),ind_loop(4)
        do i=ind_loop(1),ind_loop(2)

          l = indmtr(i  ,j  ,k  )
          l0= indcg(i  ,j  ,k  ) ! x(i  , j  , k  )
          l1= l0 + incj               ! x(i  , j+1, k  )
          l2= l0 + inck               ! x(i  , j  , k+1)
          l3= l1 + inck               ! x(i  , j+1, k+1)
          l4= l0 + inci               ! x(i+1, j  , k  )
          l5= l1 + inci               ! x(i+1, j+1, k  )
          l6= l2 + inci               ! x(i+1, j  , k+1)
          l7= l3 + inci               ! x(i+1, j+1, k+1)

          lip = l + inci_mtr              
          ljp = l + incj_mtr
          lkp = l + inck_mtr

          ax = .25*( ( x(l4) + x(l5) + x(l7) +x(l6))*ti(lip,1)
     &              -( x(l0) + x(l1) + x(l3) +x(l2))*ti(l  ,1)
     &              +( x(l1) + x(l3) + x(l5) +x(l7))*tj(ljp,1)
     &              -( x(l0) + x(l2) + x(l4) +x(l6))*tj(l  ,1)
     &              +( x(l2) + x(l3) + x(l6) +x(l7))*tk(lkp,1)
     &              -( x(l0) + x(l1) + x(l4) +x(l5))*tk(l  ,1))

          ay = .25*( ( y(l4) + y(l5) + y(l7) +y(l6))*ti(lip,2)
     &              -( y(l0) + y(l1) + y(l3) +y(l2))*ti(l  ,2)
     &              +( y(l1) + y(l3) + y(l5) +y(l7))*tj(ljp,2)
     &              -( y(l0) + y(l2) + y(l4) +y(l6))*tj(l  ,2)
     &              +( y(l2) + y(l3) + y(l6) +y(l7))*tk(lkp,2)
     &              -( y(l0) + y(l1) + y(l4) +y(l5))*tk(l  ,2))

          az = .25*( ( z(l4) + z(l5) + z(l7) +z(l6))*ti(lip,3)
     &              -( z(l0) + z(l1) + z(l3) +z(l2))*ti(l  ,3)
     &              +( z(l1) + z(l3) + z(l5) +z(l7))*tj(ljp,3)
     &              -( z(l0) + z(l2) + z(l4) +z(l6))*tj(l  ,3)
     &              +( z(l2) + z(l3) + z(l6) +z(l7))*tk(lkp,3)
     &              -( z(l0) + z(l1) + z(l4) +z(l5))*tk(l  ,3))

          vol(l) = (ax+ay+az)/3.
        enddo
        enddo
        enddo

      ELSEIF(param_int( ITYPZONE ).eq.1) THEN !Domaine 3D avec une direction homogene k: traitement facette i et j

        !calcul volume
        do k=ind_loop(5),ind_loop(6)
        do j=ind_loop(3),ind_loop(4)
        do i=ind_loop(1),ind_loop(2)

          l = indmtr(i  ,j  ,k  )
          l0= indcg(i  ,j  ,k  ) ! x(i  , j  , k  )
          l1= l0 + incj               ! x(i  , j+1, k  )
          l2= l0 + inck               ! x(i  , j  , k+1)
          l3= l1 + inck               ! x(i  , j+1, k+1)
          l4= l0 + inci               ! x(i+1, j  , k  )
          l5= l1 + inci               ! x(i+1, j+1, k  )
          l6= l2 + inci               ! x(i+1, j  , k+1)
          l7= l3 + inci               ! x(i+1, j+1, k+1)

          lip = l + inci_mtr              
          ljp = l + incj_mtr
          lkp = l + inck_mtr

          ax = .25*( ( x(l4) + x(l5) + x(l7) +x(l6))*ti(lip,1)
     &              -( x(l0) + x(l1) + x(l3) +x(l2))*ti(l  ,1)
     &              +( x(l1) + x(l3) + x(l5) +x(l7))*tj(ljp,1)
     &              -( x(l0) + x(l2) + x(l4) +x(l6))*tj(l  ,1))

          ay = .25*( ( y(l4) + y(l5) + y(l7) +y(l6))*ti(lip,2)
     &              -( y(l0) + y(l1) + y(l3) +y(l2))*ti(l  ,2)
     &              +( y(l1) + y(l3) + y(l5) +y(l7))*tj(ljp,2)
     &              -( y(l0) + y(l2) + y(l4) +y(l6))*tj(l  ,2))

          az = .25*( ( z(l2) + z(l3) + z(l6) +z(l7))*tk(lkp,1)
     &              -( z(l0) + z(l1) + z(l4) +z(l5))*tk(l  ,1))

          vol(l) = (ax+ay+az)/3.
        enddo
        enddo
        enddo

      ELSEIF(param_int( ITYPZONE ).eq.2) THEN !Domaine 3D cartesien

        !calcul volume
        do k=ind_loop(5),ind_loop(6)
        do j=ind_loop(3),ind_loop(4)
        do i=ind_loop(1),ind_loop(2)

          l = indmtr(i  ,j  ,k  )
          l0= indcg(i  ,j  ,k  ) ! x(i  , j  , k  )
          l1= l0 + incj               ! x(i  , j+1, k  )
          l2= l0 + inck               ! x(i  , j  , k+1)
          l3= l1 + inck               ! x(i  , j+1, k+1)
          l4= l0 + inci               ! x(i+1, j  , k  )
          l5= l1 + inci               ! x(i+1, j+1, k  )
          l6= l2 + inci               ! x(i+1, j  , k+1)
          l7= l3 + inci               ! x(i+1, j+1, k+1)

          lip = l + inci_mtr              
          ljp = l + incj_mtr
          lkp = l + inck_mtr

          ax = .25*( ( x(l4) + x(l5) + x(l7) +x(l6))*ti(lip,1)
     &              -( x(l0) + x(l1) + x(l3) +x(l2))*ti(l  ,1))

          ay = .25*( ( y(l1) + y(l3) + y(l5) +y(l7))*tj(ljp,1)
     &              -( y(l0) + y(l2) + y(l4) +y(l6))*tj(l  ,1))

          az = .25*( ( z(l2) + z(l3) + z(l6) +z(l7))*tk(lkp,1)
     &              -( z(l0) + z(l1) + z(l4) +z(l5))*tk(l  ,1))

          vol(l) = (ax+ay+az)/3.
        enddo
        enddo
        enddo
      ELSE !Domaine 2D: neq_k=0

        !calcul volume
        do k=ind_loop(5),ind_loop(6)
        do j=ind_loop(3),ind_loop(4)
        do i=ind_loop(1),ind_loop(2)

          l = indmtr(i  ,j  ,k  )
          l0= indcg(i  ,j  ,k  ) ! x(i  , j  , k  )
          l1= l0 + incj               ! x(i  , j+1, k  )
          l2= l0 + inck               ! x(i  , j  , k+1)
          l3= l1 + inck               ! x(i  , j+1, k+1)
          l4= l0 + inci               ! x(i+1, j  , k  )
          l5= l1 + inci               ! x(i+1, j+1, k  )
          l6= l2 + inci               ! x(i+1, j  , k+1)
          l7= l3 + inci               ! x(i+1, j+1, k+1)

          lip = l + inci_mtr              
          ljp = l + incj_mtr
          lkp = l + inck_mtr

          ax = .25*( ( x(l4) + x(l5) + x(l7) +x(l6))*ti(lip,1)
     &              -( x(l0) + x(l1) + x(l3) +x(l2))*ti(l  ,1)
     &              +( x(l1) + x(l3) + x(l5) +x(l7))*tj(ljp,1)
     &              -( x(l0) + x(l2) + x(l4) +x(l6))*tj(l  ,1))

          ay = .25*( ( y(l4) + y(l5) + y(l7) +y(l6))*ti(lip,2)
     &              -( y(l0) + y(l1) + y(l3) +y(l2))*ti(l  ,2)
     &              +( y(l1) + y(l3) + y(l5) +y(l7))*tj(ljp,2)
     &              -( y(l0) + y(l2) + y(l4) +y(l6))*tj(l  ,2))

          tkz10=(x(l4)-x(l0))*(y(l5)-y(l4))-(y(l4)-y(l0))*(x(l5)-x(l4))
          tkz11=(x(l1)-x(l5))*(y(l0)-y(l1))-(y(l1)-y(l5))*(x(l0)-x(l1))
          tk1  = .5*(tkz10+tkz11)  !=normal*surf (k)
 
          az = .25*(z(l2)+z(l3)+z(l6)+z(l7)-z(l0)-z(l1)-z(l4)-z(l5))*tk1

          vol(l) = (ax+ay+az)/3.
        enddo
        enddo
        enddo

      ENDIF!neq_k

      end
