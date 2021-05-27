c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_wallexchange(idir,lrhs, neq_mtr,
     &                          mobile_coef,c4,c5,c6,
     &                          param_int,param_real, ind_loop,
     &                          ventijk, tijk, rop, xmut)
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

      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E xmut   (param_int(NDIMDX     ), 2                   )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E mobile_coef,c4,c5,c6, param_real(0:*)
C Var local
      INTEGER_E  inc2,inc3,li1,li2,l,iref,jref,kref,lij,lr,lp,
     & njf,nkf,ldjr,ic,jc,kc,i,j,k,li3,ldp,kc_vent,sample,
     & incijk,ldjr2,lsample,lt,m,lmtr,iter,l0,ldjr3,ldjr4,
     & l1,l2,incj,inck

      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & ventx,venty,ventz,r,u,v,w,ut,vt,wt,ua,va,wa,s_1,qn,
     & nx, ny,nz,utau2,seuil,unorm1,qen,u_int,c0,c1,c2,c3,
     & dist,utau,unorm,surf,aa,yplus,g1,g2,g3,f,tp,fp,tauw,
     & Twall,uplus,pr,prtur,ratio,dt,pf,Tplus,qwall,cp,cond,
     & wall_ut,wall_vt,wall_wt,wall_int,co,si,rot1,rot2,ra,un,vn,wn,umod

#include "FastS/formule_param.h"
#include "FastS/formule_vent_param.h"
#include "FastS/formule_mtr_param.h"



c      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop vis=', ind_loop

      sample = param_int(WM_SAMPLING)


      c0   = 1./c6
      c1   =-(c4 + c5)*c0
      c2   =- c6*c0
      c3   = (2.- c5- c4)*c0

c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      c_ale=c_ale*2.*mobile_coef
      if(lrhs.eq.1) c_ale = 0.

      IF (idir.eq.1) THEN

       if(param_int(NEQ).eq.5) then

          do 100 k = ind_loop(5), ind_loop(6)
          do 100 j = ind_loop(3), ind_loop(4)

             i = ind_loop(2)

             l      = inddm(  i        , j,  k ) 
             ldjr   = inddm(  i+1      , j,  k )
             ldjr2  = inddm(  i+2      , j,  k )
             lsample= inddm(  i+sample , j,  k )
             ldp    = indven( i+1      , j,  k )
             lmtr   = indmtr( i+1      , j,  k )

#include     "FastS/BC/BCWallExchange.for"

             l0   = l
             l1   = l + 1
             l2   = l + 2
             do i = ind_loop(1), ind_loop(2)-1

                  l    = inddm( i , j,  k ) 
#include          "FastS/BC/BC_nextrank.for"
             enddo

100    continue

       else
        ! a faire
       endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

          do 120 k = ind_loop(5), ind_loop(6)
          do 120 j = ind_loop(3), ind_loop(4)

              i      = ind_loop(1)
              l      = inddm(  i         , j, k ) 
              ldjr   = inddm(  i -1      , j, k )
              ldjr2  = inddm(  i -2      , j, k )
              lsample= inddm(  i -sample , j, k )
              ldp    = indven( i         , j, k )
              lmtr   = indmtr( i         , j, k )

#include     "FastS/BC/BCWallExchange.for"

               l0   = l
               l1   = l-1
               l2   = l-2
               do i = ind_loop(1)+1, ind_loop(2)

                 l    = inddm(  i , j, k ) 
#include        "FastS/BC/BC_nextrank.for"
               enddo

120    continue


      ELSEIF (idir.eq.3) THEN

          incj = param_int(NIJK)

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
                !do j = ind_loop(4)+1, ind_loop(4)+2 
                do j = ind_loop(4)+1, ind_loop(4)+3 
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  !j    = ind_loop(4)
                  l    = inddm( i ,  j              , k )
                  ldjr = inddm( i ,  ind_loop(4) +1 , k )
                  !ldjr3= inddm( i ,  ind_loop(4) +3 , k )
                  ldjr3= inddm( i ,  ind_loop(4) +4 , k )
                  ldjr4= inddm( i ,  ind_loop(4) +5 , k )
                  lmtr = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp  = indven(i ,  ind_loop(4) +1 , k )

                  !ldjr   = l +   incj
                  !ldjr2  = l + 2*incj
                  lsample= ldjr + (sample-1)*incj
#include          "FastS/BC/BCWallExchange.for"
                enddo !i
                enddo !j

                do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2) 

                     l    = inddm( i ,    j            , k )
                     l0   = inddm( i ,  ind_loop(4)+1  , k )

                      rop(l,1) = rop(l0,1) 
                      rop(l,2) = rop(l0,2)
                      rop(l,3) = rop(l0,3)
                      rop(l,4) = rop(l0,4) 
                      rop(l,5) = rop(l0,5)
c                     l    = inddm( i ,    j            , k )
c                     l0   = inddm( i ,  ind_loop(4)    , k )
c                     l1   = l0 +   incj
c                     l2   = l0 + 2*incj
CC#include             "FastS/BC/BC_nextrank.for"
                  enddo!i
                enddo !j
             enddo !k

       else
       endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

          incj = -param_int(NIJK)

             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  j      = ind_loop(3)
                  l      = inddm( i ,  j  , k )
                  lmtr   = indmtr(i ,  j  , k )
                  ldp    = indven(i ,  j  , k )
                  ldjr   = l +   incj
                  ldjr2  = l + 2*incj
                  lsample= l + sample*incj

#include          "FastS/BC/BCWallExchange.for"
                enddo !i

                do  j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2) 

                     l    = inddm( i ,    j          , k )
                     l0   = inddm( i ,  ind_loop(3)  , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_nextrank.for"
                  enddo!i
                enddo !j
             enddo !k

      ELSEIF (idir.eq.5) THEN

          inck = param_int(NIJK)*param_int(NIJK+1)

             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)

                  k    = ind_loop(6)
                  l    = inddm( i , j,  k   )
                  lmtr = indmtr(i , j,  k+1 )
                  ldp  = indven(i , j,  k+1 )
                  ldjr   = l +   inck
                  ldjr2  = l + 2*inck
                  lsample= l + sample*inck
#include        "FastS/BC/BCWallExchange.for"
               enddo
             enddo

             do  k = ind_loop(5), ind_loop(6)-1
                do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
                   do  i = ind_loop(1), ind_loop(2)

                     l    = inddm( i , j,  k          )
                     l0   = inddm( i , j, ind_loop(6) )
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank.for"
                   enddo
               enddo
             enddo

      ELSE 

          inck = -param_int(NIJK)*param_int(NIJK+1)

             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)

                  k      = ind_loop(5)
                  l      = inddm( i , j,  k   )
                  lmtr   = indmtr(i , j,  k   )
                  ldp    = indven(i , j,  k   )
                  ldjr   = l +   inck
                  ldjr2  = l + 2*inck
                  lsample= l + sample*inck
#include        "FastS/BC/BCWallExchange.for"
               enddo
             enddo

             do  k = ind_loop(5)+1, ind_loop(6)
                do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
                   do  i = ind_loop(1), ind_loop(2)

                     l    = inddm( i , j,  k          )
                     l0   = inddm( i , j, ind_loop(5) )
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank.for"
                   enddo
               enddo
             enddo

      ENDIF !idir


      END
