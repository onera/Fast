c***********************************************************************
c     $Date: 2011-06-20 11:10:00 +0200 (lun 20 jun 2011) $
c     $Revision: 58 $
c     $Author: NicolasAlferez $
c***********************************************************************
      subroutine move_domx(ndom, param_int, param_real,
     &                      ind_loop,
     &                      x,y,z,xinit,yinit,zinit)
c***********************************************************************
c_P                          O N E R A
c
c_DC  DATE_C : Janvier 2011  --  ALFEREZ
c
c_DM  DATE_M : 17 Jan, 13:57
c
c_U   USER : ALFEREZ
c
c     ACT
c_A   Appel de la routine de calcul du torseur du mouvement relatif
c_A   du maillage
c     INP
c_I    temps  :  temps a l'iteration n (pas de temps)
c_I              ou n+1/2 
c     OUT : xm,ym, zm: coordonnee maillage a l'instant t pour extraction data ou flow
c     I/O
c
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), param_int(0:*)

      REAL_E x( param_int(NDIMDX_XYZ) ),y( param_int(NDIMDX_XYZ) ),
     &       z( param_int(NDIMDX_XYZ) )
      REAL_E xinit(param_int(NDIMDX_XYZ)),yinit( param_int(NDIMDX_XYZ)),
     &       zinit(param_int(NDIMDX_XYZ)), param_real(0:*)

C var loc  
      INTEGER_E l,inck,k,i,j,incmax,i2save,j2save,k2save
      INTEGER_E translation_pur
      REAL_E rot(4,3), vtrans(3)

#include "FastS/formule_xyz_param.h"

      if(ind_loop(1).gt.ind_loop(2)) goto 10
      if(ind_loop(3).gt.ind_loop(4)) goto 10
      if(ind_loop(5).gt.ind_loop(6)) goto 10

       vtrans=0.
       !om modifie la position en fonction de la loi horaire
       call move( rot, param_real, translation_pur) 

#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
      do l = 1,ndimt_xyz
#else
      do k=ind_loop(5),ind_loop(6)
        do j=ind_loop(3),ind_loop(4)
          do i=ind_loop(1),ind_loop(2)

          l= indcg(i  ,j  ,k     ) ! (i  , j  , k  )
#endif
      x(l) = param_real(ROT_CENTER  )
     &                    + rot(1,1)*(xinit(l)-param_real(ROT_CENTER  ))
     &                    + rot(1,2)*(yinit(l)-param_real(ROT_CENTER+1))
     &                    + rot(1,3)*(zinit(l)-param_real(ROT_CENTER+2))

      y(l) = param_real(ROT_CENTER+1) 
     &                    + rot(2,1)*(xinit(l)-param_real(ROT_CENTER  ))
     &                    + rot(2,2)*(yinit(l)-param_real(ROT_CENTER+1))
     &                    + rot(2,3)*(zinit(l)-param_real(ROT_CENTER+2))

      z(l) = param_real(ROT_CENTER+2)
     &                    + rot(3,1)*(xinit(l)-param_real(ROT_CENTER  ))
     &                    + rot(3,2)*(yinit(l)-param_real(ROT_CENTER+1))
     &                    + rot(3,3)*(zinit(l)-param_real(ROT_CENTER+2))
 
#ifndef E_SCALAR_COMPUTER
      enddo
#else
        enddo
        enddo
        enddo
#endif

 10   continue

      end
