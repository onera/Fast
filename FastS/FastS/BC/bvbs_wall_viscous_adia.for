c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_wall_viscous_adia(idir,lrhs, neq_mtr, mobile_coef,
     &                                  param_int, ind_loop,
     &                                  ventijk, tijk, rop)
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
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E mobile_coef
C Var local
      INTEGER_E  inc2,inc3,li1,li2,l,iref,jref,kref,lij,lr,lp,
     & njf,nkf,ldjr,ic,jc,kc,i,j,k,li3,ldp,kc_vent
      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & ventx,venty,ventz,r,u,v,w,ut,vt,wt,ua,va,wa,s_1,qn

#include "FastS/formule_param.h"
#include "FastS/formule_vent_param.h"



c      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop vis=', ind_loop


c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      c_ale=c_ale*2.*mobile_coef
      if(lrhs.eq.1) c_ale = 0.

      IF (idir.eq.1) THEN

       iref = 2*ind_loop(2) + 1


       if(param_int(NEQ).eq.5) then

          do 100 k = ind_loop(5), ind_loop(6)
          do 100 j = ind_loop(3), ind_loop(4)

            do 100 i = ind_loop(1), ind_loop(2)

             l    = inddm( i              , j,  k ) 
             ldjr = inddm(  iref - i      , j,  k )
             ldp  = indven( ind_loop(2)+1 , j,  k )
#include     "FastS/BC/BCWallViscous.for"
100    continue

       else

          do 110 k = ind_loop(5), ind_loop(6)
          do 110 j = ind_loop(3), ind_loop(4)
          do 110 i = ind_loop(1), ind_loop(2)

             l    = inddm( i              , j,  k ) 
             ldjr = inddm(  iref - i      , j,  k )
             ldp  = indven( ind_loop(2)+1 , j,  k )

#include     "FastS/BC/BCWallViscous_SA.for"
110       continue
        endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

       iref = 2*ind_loop(1) - 1

       if(param_int(NEQ).eq.5) then


          do 120 k = ind_loop(5), ind_loop(6)
          do 120 j = ind_loop(3), ind_loop(4)
          do 120 i = ind_loop(1), ind_loop(2)

              l    = inddm(  i           , j, k ) 
              ldjr = inddm(  iref - i    , j, k )
              ldp  = indven( ind_loop(1) , j, k )
#include     "FastS/BC/BCWallViscous.for"
120    continue
       else

          do 130 k = ind_loop(5), ind_loop(6)
          do 130 j = ind_loop(3), ind_loop(4)
          do 130 i = ind_loop(1), ind_loop(2)

              l    = inddm(  i           , j, k ) 
              ldjr = inddm(  iref - i    , j, k )
              ldp  = indven( ind_loop(1) , j, k )
#include     "FastS/BC/BCWallViscous_SA.for"
130       continue
        endif !param_int(NEQ)


      ELSEIF (idir.eq.3) THEN

       jref = 2*ind_loop(4) + 1

       if(param_int(NEQ).eq.5) then

          do 200 k = ind_loop(5), ind_loop(6)
          do 200 j = ind_loop(3), ind_loop(4)

            lij =       inddm( ind_loop(1) , j             , k )
            lr  = lij - inddm( ind_loop(1) , jref - j      , k )
            lp  = lij - indven(ind_loop(1) , ind_loop(4)+1 , k )

!DEC$ IVDEP
            do 200 l = lij, lij + ind_loop(2) - ind_loop(1)

                ldjr = l - lr
                ldp  = l - lp
#include     "FastS/BC/BCWallViscous.for"
200    continue
       else

          do 210 k = ind_loop(5), ind_loop(6)
          do 210 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j             , k )
            lr  = lij - inddm( ind_loop(1) , jref - j      , k )
            lp  = lij - indven(ind_loop(1) , ind_loop(4)+1 , k )

!DEC$ IVDEP
            do 210 l = lij, lij + ind_loop(2) - ind_loop(1)

                ldjr = l - lr
                ldp  = l - lp
#include     "FastS/BC/BCWallViscous_SA.for"
210       continue
        endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

       jref = 2*ind_loop(3) - 1

       if(param_int(NEQ).eq.5) then

          do 220 k = ind_loop(5), ind_loop(6)
          do 220 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j             , k )
            lr  = lij - inddm( ind_loop(1) , jref - j      , k )
            lp  = lij - indven(ind_loop(1) , ind_loop(3)   , k )
!DEC$ IVDEP
            do 220 l = lij, lij + ind_loop(2) - ind_loop(1)
              ldjr = l - lr
              ldp  = l - lp
#include     "FastS/BC/BCWallViscous.for"
220    continue

       else

          do 230 k = ind_loop(5), ind_loop(6)
          do 230 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j           , k )
            lr  = lij - inddm( ind_loop(1) , jref - j    , k )
            lp  = lij - indven(ind_loop(1) , ind_loop(3) , k )
!DEC$ IVDEP
            do 230 l = lij, lij + ind_loop(2) - ind_loop(1)
              ldjr = l - lr
              ldp  = l - lp
#include     "FastS/BC/BCWallViscous_SA.for"
230       continue
        endif !param_int(NEQ)


      ELSEIF (idir.eq.5) THEN

       kref = 2*ind_loop(6) + 1

       if(param_int(NEQ).eq.5) then

          do 300 k = ind_loop(5), ind_loop(6)
          do 300 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j  , k             )
            lr  = lij - inddm( ind_loop(1) , j  , kref - k      )
            lp  = lij - indven(ind_loop(1) , j ,  ind_loop(6)+1 )
!DEC$ IVDEP
            do 300 l = lij, lij + ind_loop(2) - ind_loop(1)

              ldjr = l - lr
              ldp  = l - lp
#include     "FastS/BC/BCWallViscous.for"
300    continue
       else
          do 310 k = ind_loop(5), ind_loop(6)
          do 310 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j  , k             )
            lr  = lij - inddm( ind_loop(1) , j  , kref - k      )
            lp  = lij - indven(ind_loop(1) , j ,  ind_loop(6)+1 )
!DEC$ IVDEP
            do 310 l = lij, lij + ind_loop(2) - ind_loop(1)

              ldjr = l - lr
              ldp  = l - lp
#include     "FastS/BC/BCWallViscous_SA.for"
310       continue
        endif !param_int(NEQ)


      ELSE 

       kref = 2*ind_loop(5) - 1

       if(param_int(NEQ).eq.5) then


          do 320 k = ind_loop(5), ind_loop(6)
          do 320 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j  , k           )
            lr  = lij - inddm( ind_loop(1) , j  , kref - k    )
            lp  = lij - indven(ind_loop(1) , j ,  ind_loop(5) )
!DEC$ IVDEP
            do 320 l = lij, lij + ind_loop(2) - ind_loop(1)

              ldjr = l - lr
              ldp  = l - lp
#include     "FastS/BC/BCWallViscous.for"
320    continue
       else

          do 330 k = ind_loop(5), ind_loop(6)
          do 330 j = ind_loop(3), ind_loop(4)
            lij =       inddm( ind_loop(1) , j  , k           )
            lr  = lij - inddm( ind_loop(1) , j  , kref - k    )
            lp  = lij - indven(ind_loop(1) , j ,  ind_loop(5) )
!DEC$ IVDEP
            do 330 l = lij, lij + ind_loop(2) - ind_loop(1)

              ldjr = l - lr
              ldp  = l - lp
#include     "FastS/BC/BCWallViscous_SA.for"
330       continue
        endif !param_int(NEQ)

      ENDIF !idir


      END
