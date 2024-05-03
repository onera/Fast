c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_wallexchange(idir,lrhs, neq_mtr,
     &                          mobile_coef,c4,c5,c6,
     &                          param_int,param_real, ind_loop,
     &                          sizetab, inc_bc,
     &                          x, y,z, ventijk, tijk, rop, xmut,
     &                          moy, param_wmles)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c     Paroi glissante
c     VAL
c_V    Optimisation NEC
c
c     COM
c_C    MODIFIER bvas3.f (passage de ird1,ird2a,ird2b ?????)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E idir,lrhs, neq_mtr, ind_loop(6), param_int(0:*)
      INTEGER_E sizetab, Nechant, inc_bc(3)

      REAL_E x( param_int(NDIMDX_XYZ) ),y( param_int(NDIMDX_XYZ) ),
     &       z( param_int(NDIMDX_XYZ) )

      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E xmut   (param_int(NDIMDX     ), 2                   )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E moy    (sizetab, 5 )
      REAL_E mobile_coef,c4,c5,c6, param_real(0:*),param_wmles(*)
C Var local
      INTEGER_E  inc2,inc3,li1,li2,l,iref,jref,kref,lij,lr,lp,
     & njf,nkf,ldjr,ic,jc,kc,i,j,k,li3,ldp,kc_vent,sample,
     & lsample,lt,m,lmtr,iter,l0,l1,l2,incj,inck,lghost

      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & ventx,venty,ventz,r,u,v,w,ut,vt,wt,ua,va,wa,s_1,qn,dx,dy,dz,
     & nx, ny,nz,utau2,seuil,unorm1,qen,u_int,c0,c1,c2,c3,
     & dist,utau,unorm,surf,aa,yplus,g1,g2,g3,f,tp,fp,tauw,
     & Twall,uplus,pr,prtur,ratio,dt,pf,Tplus,qwall,cp,cond,vmoy_bc,
     & wall_ut,wall_vt,wall_wt,wall_int,co,si,rot1,rot2,ra,un,vn,wn,
     & umod,rescale, gam,cv,prandlt,pwall,tol,cv1,
     & expy,nutcible,nu,a,b,p,q,ap,bp,dp,delta,delta0,delta1,
     & y1,y2,z1,p1,p2,kappa,superdelta,root,b0,cmus1,temp01,coesut

#include "FastS/formule_param.h"
#include "FastS/formule_xyz_param.h"
#include "FastS/formule_vent_param.h"
#include "FastS/formule_mtr_param.h"

      !-----Variables physiques
      gam     = param_real( GAMMA )
      cv      = param_real( CVINF )
      prandlt = param_real( PRANDT )

      cmus1  =    param_real( CS )
      temp01 = 1./param_real( TEMP0 )
      coesut =    param_real( XMUL0 ) * (1.+cmus1*temp01)*sqrt(temp01)


      tol   = 1.e-12
      cv1   = 7.1
      kappa = 0.4

      incj = param_int(NIJK)
      inck = param_int(NIJK)* param_int(NIJK+1)

c      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop vis=', ind_loop, inc_bc

      sample = param_int(WM_SAMPLING)

c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      c_ale=c_ale*2.*mobile_coef
      if(lrhs.eq.1) c_ale = 0.

      IF (idir.eq.1) THEN

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

             !1ere rangee reelle au dessus paroi
             i      = ind_loop(2)+1
             l      = inddm( i              , j, k )
             ldjr   = inddm( ind_loop(2) +1 , j, k )
             lmtr   = indmtr(ind_loop(2) +1 , j, k )
             ldp    = indven(ind_loop(2) +1 , j, k )
             lsample= ldjr + (sample-1)
#include     "FastS/BC/BCWallExchange.for"

             !extrap 1ere rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(2)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)

             !2eme rangee reelle au dessus paroi
             i      = ind_loop(2)+2
             l      = inddm( i              , j, k )
#include     "FastS/BC/BCWallExchange.for"
             !extrap 2eme rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(1)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)
          enddo
          enddo
       else
         do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

             !1ere rangee reelle au dessus paroi
             i      = ind_loop(2)+1
             l      = inddm( i              , j, k )
             ldjr   = inddm( ind_loop(2) +1 , j, k )
             lmtr   = indmtr(ind_loop(2) +1 , j, k )
             ldp    = indven(ind_loop(2) +1 , j, k )
             lsample= ldjr + (sample-1)
#include     "FastS/BC/BCWallExchangeSA.for"

             !extrap 1ere rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(2)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)
             rop(lghost,6) = rop(ldjr,6)

             !2eme rangee reelle au dessus paroi
             i      = ind_loop(2)+2
             l      = inddm( i              , j, k )
#include     "FastS/BC/BCWallExchangeSA.for"
             !extrap 2eme rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(1)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)
             rop(lghost,6) = rop(ldjr,6)
          enddo
          enddo
       endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
             !1ere rangee reelle au dessus paroi
             i      = ind_loop(1)-1
             l      = inddm( i              , j, k )
             ldjr   = inddm( ind_loop(1) -1 , j, k )
             lmtr   = indmtr(ind_loop(1) -1 , j, k )
             ldp    = indven(ind_loop(1) -1 , j, k )
             lsample= ldjr - (sample-1)
#include     "FastS/BC/BCWallExchange.for"

             !extrap 1ere rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(1)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)

             !2eme rangee reelle au dessus paroi
             i      = ind_loop(1)-2
             l      = inddm( i              , j, k )
#include     "FastS/BC/BCWallExchange.for"
             !extrap 2eme rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(2)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)
          enddo
          enddo
       else
          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
             !1ere rangee reelle au dessus paroi
             i      = ind_loop(1)-1
             l      = inddm( i              , j, k )
             ldjr   = inddm( ind_loop(1) -1 , j, k )
             lmtr   = indmtr(ind_loop(1) -1 , j, k )
             ldp    = indven(ind_loop(1) -1 , j, k )
             lsample= ldjr - (sample-1)
#include     "FastS/BC/BCWallExchangeSA.for"

             !extrap 1ere rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(1)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)
             rop(lghost,6) = rop(ldjr,6)

             !2eme rangee reelle au dessus paroi
             i      = ind_loop(1)-2
             l      = inddm( i              , j, k )
#include     "FastS/BC/BCWallExchangeSA.for"
             !extrap 2eme rangee ghost au dessous paroi (pour visu)
             i      = ind_loop(2)
             lghost = inddm( i              , j, k )
             rop(lghost,1) = rop(ldjr,1) 
             rop(lghost,2) = rop(ldjr,2)
             rop(lghost,3) = rop(ldjr,3)
             rop(lghost,4) = rop(ldjr,4) 
             rop(lghost,5) = rop(ldjr,5)
             rop(lghost,6) = rop(ldjr,6)
          enddo
          enddo
       endif

      ELSEIF (idir.eq.3) THEN

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
                do j = ind_loop(4)+1, ind_loop(4)+2 
                !do j = ind_loop(4)+1, ind_loop(4)+3 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j              , k )
                  ldjr   = inddm( i ,  ind_loop(4) +1 , k )
                  lmtr   = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp    = indven(i ,  ind_loop(4) +1 , k )
                  lsample= ldjr + (sample-1)*incj
#include          "FastS/BC/BCWallExchange.for"

                  !remplissage ghost pour visu
                  lghost   = inddm( i , j-2  , k )
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                enddo !i
                enddo !j
             enddo !k

       else
             do  k = ind_loop(5), ind_loop(6)
                do j = ind_loop(4)+1, ind_loop(4)+2 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j              , k )
                  ldjr   = inddm( i ,  ind_loop(4) +1 , k )
                  lmtr   = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp    = indven(i ,  ind_loop(4) +1 , k )
                  lsample= ldjr + (sample-1)*incj

#include          "FastS/BC/BCWallExchangeSA.for"
                  !remplissage ghost pour visu
                  lghost   = inddm( i , j-2  , k )
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                  rop(lghost,6) = rop(ldjr,6)
                enddo !i
                enddo !j
             enddo !k
       endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
                do j = ind_loop(3)-2, ind_loop(3)-1
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j              , k )
                  ldjr   = inddm( i ,  ind_loop(3) -1 , k )
                  lmtr   = indmtr(i ,  ind_loop(3) -1 , k )
                  ldp    = indven(i ,  ind_loop(3) -1 , k )
                  lsample= ldjr - (sample-1)*incj
#include          "FastS/BC/BCWallExchange.for"

                  !remplissage ghost pour visu
                  lghost   = inddm( i , j+2  , k )
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                enddo !i
                enddo !j
             enddo !k

       else
             do  k = ind_loop(5), ind_loop(6)
                do j = ind_loop(3)-2, ind_loop(3)-1
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j              , k )
                  ldjr   = inddm( i ,  ind_loop(3) -1 , k )
                  lmtr   = indmtr(i ,  ind_loop(3) -1 , k )
                  ldp    = indven(i ,  ind_loop(3) -1 , k )
                  lsample= ldjr - (sample-1)*incj

#include          "FastS/BC/BCWallExchangeSA.for"
                  !remplissage ghost pour visu
                  lghost   = inddm( i , j+2  , k )
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                  rop(lghost,6) = rop(ldjr,6)
                enddo !i
                enddo !j
             enddo !k
       endif !param_int(NEQ)

      ELSEIF (idir.eq.5) THEN

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(6)+1, ind_loop(6)+2
                do j = ind_loop(3), ind_loop(4) 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j,      k         )
                  ldjr   = inddm( i ,  j, ind_loop(6) +1 )
                  lmtr   = indmtr(i ,  j, ind_loop(6) +1 )
                  ldp    = indven(i ,  j, ind_loop(6) +1 )
                  lsample= ldjr + (sample-1)*inck
#include          "FastS/BC/BCWallExchange.for"

                  !remplissage ghost pour visu
                  lghost   = inddm( i , j,  k-2)
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                enddo !i
                enddo !j
             enddo !k

       else
             do  k = ind_loop(6)+1, ind_loop(6)+2
                do j = ind_loop(3), ind_loop(4) 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j,        k       )
                  ldjr   = inddm( i ,  j, ind_loop(6) +1 )
                  lmtr   = indmtr(i ,  j, ind_loop(6) +1 )
                  ldp    = indven(i ,  j, ind_loop(6) +1 )
                  lsample= ldjr + (sample-1)*inck

#include          "FastS/BC/BCWallExchangeSA.for"
                  !remplissage ghost pour visu
                  lghost   = inddm( i , j,  k-2)
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                  rop(lghost,6) = rop(ldjr,6)
                enddo !i
                enddo !j
             enddo !k
       endif !param_int(NEQ)

      ELSE 

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5)-2, ind_loop(5)-1
                do j = ind_loop(3), ind_loop(4) 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j,      k         )
                  ldjr   = inddm( i ,  j, ind_loop(5) -1 )
                  lmtr   = indmtr(i ,  j, ind_loop(5) -1 )
                  ldp    = indven(i ,  j, ind_loop(5) -1 )
                  lsample= ldjr - (sample-1)*inck
#include          "FastS/BC/BCWallExchange.for"

                  !remplissage ghost pour visu
                  lghost   = inddm( i , j,  k+2)
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                enddo !i
                enddo !j
             enddo !k

       else
             do  k = ind_loop(5)-2, ind_loop(5)-1
                do j = ind_loop(3), ind_loop(4) 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l      = inddm( i ,  j,        k       )
                  ldjr   = inddm( i ,  j, ind_loop(5) -1 )
                  lmtr   = indmtr(i ,  j, ind_loop(5) -1 )
                  ldp    = indven(i ,  j, ind_loop(5) -1 )
                  lsample= ldjr - (sample-1)*inck

#include          "FastS/BC/BCWallExchangeSA.for"
                  !remplissage ghost pour visu
                  lghost   = inddm( i , j,  k+2)
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                  rop(lghost,6) = rop(ldjr,6)
                enddo !i
                enddo !j
             enddo !k
       endif !param_int(NEQ)

      ENDIF !idir


      END
