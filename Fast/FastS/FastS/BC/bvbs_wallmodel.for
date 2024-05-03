c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_wallmodel(idir,lrhs, neq_mtr, mobile_coef,
     &                             c4,c5,c6,
     &                             param_int, ind_loop,
     &                             ventijk, tijk, rop, xmut)
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
      REAL_E mobile_coef, c4,c5,c6
C Var local
      INTEGER_E  inc2,inc3,li1,li2,l,iref,jref,kref,
     & njf,nkf,ldnm,ldl,ldnp,ic,jc,kc,i,j,k,li3,ldp,kc_vent,lmtr
      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,c7,
     & ventx,venty,ventz,r,u,v,w,ut,vt,wt,ua,va,wa,s_1,qn,surf

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

!      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
!      write(*,*)'nijk=', nijk
!      write(*,*)'loop=', ind_loop

       c7 = (c4-c5)/c6

c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      ! c_ale=0: paroi fixe, 1 sinon
      c_ale=c_ale*mobile_coef
      
      if(lrhs.eq.1) c_ale = 0.

      IF (idir.eq.1) THEN

       iref = 2*ind_loop(2) + 1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

            l    = inddm(  i             , j,  k ) 
            ldnm = inddm(  ind_loop(2)+1 , j,  k )
            ldl  = inddm(  ind_loop(2)+2 , j,  k )
            ldnp = inddm(  ind_loop(2)+3 , j,  k )


#include     "FastS/BC/BCWallModel.for"

          enddo
          enddo
          enddo

       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

            l    = inddm(  i             , j,  k ) 
            ldnm = inddm(  ind_loop(2)+1 , j,  k )
            ldl  = inddm(  ind_loop(2)+2 , j,  k )
            ldnp = inddm(  ind_loop(2)+3 , j,  k )

#include     "FastS/BC/BCWallModel.for"
              rop(l,6) =rop(ldnm,6)

          enddo
          enddo
          enddo
        endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

       iref = 2*ind_loop(1) - 1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

              l    = inddm(  i           , j, k ) 
              ldnm = inddm(  ind_loop(1)-1 , j,  k )
              ldl  = inddm(  ind_loop(1)-2 , j,  k )
              ldnp = inddm(  ind_loop(1)-3 , j,  k )

#include     "FastS/BC/BCWallModel.for"

          enddo
          enddo
          enddo
       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

              l    = inddm(  i           , j, k ) 
              ldnm = inddm(  ind_loop(1)-1 , j,  k )
              ldl  = inddm(  ind_loop(1)-2 , j,  k )
              ldnp = inddm(  ind_loop(1)-3 , j,  k )
#include     "FastS/BC/BCWallModel.for"
              rop(l,6) = rop(ldnm,6)

          enddo
          enddo
          enddo
        endif !param_int(NEQ)


      ELSEIF (idir.eq.3) THEN


       jref = 2*ind_loop(4) + 1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

            l    = inddm( i , j            ,  k ) 
            ldnm = inddm( i , ind_loop(4)+1,  k )
            ldl  = inddm( i , ind_loop(4)+2,  k )
            ldnp = inddm( i , ind_loop(4)+3,  k )

#include     "FastS/BC/BCWallModel.for"

          enddo
          enddo
          enddo
       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

            l    = inddm( i , j            ,  k ) 
            ldnm = inddm( i , ind_loop(4)+1,  k )
            ldl  = inddm( i , ind_loop(4)+2,  k )
            ldnp = inddm( i , ind_loop(4)+3,  k )

#include     "FastS/BC/BCWallModel.for"
              rop(l,6) = rop(ldnm,6)

          enddo
          enddo
          enddo
        endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

       jref = 2*ind_loop(3) - 1

       if(param_int(NEQ).eq.5) then


          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
            do i = ind_loop(1), ind_loop(2)

              l    = inddm( i , j            ,  k ) 
              ldnm = inddm( i , ind_loop(3)-1,  k )
              ldl  = inddm( i , ind_loop(3)-2,  k )
              ldnp = inddm( i , ind_loop(3)-3,  k )

#include     "FastS/BC/BCWallModel.for"

            enddo
          enddo
          enddo
       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

              l    = inddm( i , j            ,  k ) 
              ldnm = inddm( i , ind_loop(3)-1,  k )
              ldl  = inddm( i , ind_loop(3)-2,  k )
              ldnp = inddm( i , ind_loop(3)-3,  k )
#include     "FastS/BC/BCWallModel.for"
              rop(l,6) = rop(ldnm,6)

          enddo
          enddo
          enddo

        endif !param_int(NEQ)


      ELSEIF (idir.eq.5) THEN

       kref = 2*ind_loop(6) + 1

       if(param_int(NEQ).eq.5) then

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

            l    = inddm( i , j,  k            ) 
            ldnm = inddm( i , j, ind_loop(6)+1 )
            ldl  = inddm( i , j, ind_loop(6)+2 )
            ldnp = inddm( i , j, ind_loop(6)+3 )
#include     "FastS/BC/BCWallModel.for"

          enddo
          enddo
          enddo
       else
          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

            l    = inddm( i , j,  k             ) 
            ldnm = inddm( i , j, ind_loop(6)+1 )
            ldl  = inddm( i , j, ind_loop(6)+2 )
            ldnp = inddm( i , j, ind_loop(6)+3 )
#include    "FastS/BC/BCWallModel.for"
            rop(l,6) = rop(ldnm,6)

          enddo
          enddo
          enddo

        endif !param_int(NEQ)

      ELSE 

       kref = 2*ind_loop(5) - 1

       if(param_int(NEQ).eq.5) then


          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
            do i = ind_loop(1), ind_loop(2)

              l    = inddm( i , j,  k            ) 
              ldnm = inddm( i , j, ind_loop(5)-1 )
              ldl  = inddm( i , j, ind_loop(5)-2 )
              ldnp = inddm( i , j, ind_loop(5)-3 )
#include     "FastS/BC/BCWallModel.for"
            enddo
           enddo
           enddo
       else

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)
          do i = ind_loop(1), ind_loop(2)

              l    = inddm( i , j,  k            ) 
              ldnm = inddm( i , j, ind_loop(5)-1 )
              ldl  = inddm( i , j, ind_loop(5)-2 )
              ldnp = inddm( i , j, ind_loop(5)-3 )
#include     "FastS/BC/BCWallModel.for"
              rop(l,6) = rop(ldnm,6)

          enddo
          enddo
          enddo
        endif !param_int(NEQ)

      ENDIF !idir

      END
