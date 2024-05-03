c***********************************************************************
c     $Date: 2018-08-29 16:00:23 +0200 (Jeu. 29 novembre 2019) $
c     $Revision:  $
c     $Author: IvanMary  $
c***********************************************************************
      subroutine bvbs_wall_viscous_isothermal(idir,lrhs, neq_mtr,
     &                            mobile_coef,param_int, ind_loop,
     &                            param_real,
     &                            ventijk, tijk, rop,
     &                            tw,
     &                            size_data, inc_bc, size_work)
c d0x, d0y, d0z,  pa, ha, nue, -> tw
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
      INTEGER_E size_data,inc_bc(3), size_work

      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E state(param_int(NEQ))
      REAL_E param_real(0:*)
      REAL_E mobile_coef

C...  Contient les donnees CGNS pour la paroi isotherme:
      REAL_E tw(size_data)

C Var local
      INTEGER_E l,lij,lr,lp,i,j,k,l1,l2,ic,jc,kc,kc_vent,ldp,lmtr,l0,
     & incj, inck,indbci,li,ldjr,iref,jref,kref
      REAL_E ro,rho_w,u,v,w,t,nut,c6inv,c0,c1,c2,c3,roe,rue,rve,rwe,ete,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & qref1,qref2,qref3,qref4,qref5,r,p,c,ci,snorm, gamm1,
     & gamm1_1, cvinv,
     & sni,tn,qn,qen,eigv1,eigv2,eigv3,eigv4,
     & eigv5,qvar1,qvar2,qvar3,qvar4,qvar5,svar1,svar2,svar3,
     & svar4,svar5,rvar1,rvar2,rvar3,rvar4,rvar5,
     & roext,ruext,rvext,rwext,etext,roint,ruint,rvint,rwint,etint,
     & s_1,roi,rui,rvi,rwi,eti, tnx,tny, ri, roinv, sn, roe_inv,
     & mach, pref, roi0, ru,
     & ventx,venty,ventz

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

      indbci(j_1,k_1) = 1 + (j_1-inc_bc(2)) + (k_1-inc_bc(3))*inc_bc(1)

c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
c      write(*,*)'tw',tw
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      c_ale=c_ale*2.*mobile_coef
      if(lrhs.eq.1) c_ale = 0.

      IF (idir.eq.1) THEN

      iref = 2*ind_loop(2) + 1


      if(param_int(NEQ).eq.5) then

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           li   = indbci(j,  k )
           do i = ind_loop(1), ind_loop(2)

            l    = inddm( i              , j,  k )
            ldjr = inddm(  iref - i      , j,  k )
            ldp  = indven( ind_loop(2)+1 , j,  k )

#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"

           enddo
          enddo
          enddo
      else

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
            li   = indbci(j,  k )
            do i = ind_loop(1), ind_loop(2)
              l    = inddm( i              , j,  k )
              ldjr = inddm(  iref - i      , j,  k )
              ldp  = indven( ind_loop(2)+1 , j,  k )
#include     "FastS/BC/BCWallViscousIsothermal_firstrank_SA.for"

            enddo
         enddo
         enddo

       endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

      iref = 2*ind_loop(1) - 1

      if(param_int(NEQ).eq.5) then


         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
         li   = indbci(j,  k )
         do i = ind_loop(1), ind_loop(2)

             l    = inddm(  i           , j, k )
             ldjr = inddm(  iref - i    , j, k )
             ldp  = indven( ind_loop(1) , j, k )
#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"
         enddo
         enddo
         enddo

      else

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
         li   = indbci(j,  k )
         do i = ind_loop(1), ind_loop(2)

             l    = inddm(  i           , j, k )
             ldjr = inddm(  iref - i    , j, k )
             ldp  = indven( ind_loop(1) , j, k )
#include     "FastS/BC/BCWallViscousIsothermal_firstrank_SA.for"

         enddo
         enddo
         enddo
       endif !param_int(NEQ)


      ELSEIF (idir.eq.3) THEN

      jref = 2*ind_loop(4) + 1

      if(param_int(NEQ).eq.5) then

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)

           lij =       inddm( ind_loop(1) , j             , k )
           lr  = lij - inddm( ind_loop(1) , jref - j      , k )
           lp  = lij - indven(ind_loop(1) , ind_loop(4)+1 , k )

!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
               li   = indbci(l-lij+ind_loop(1),  k )
               ldjr = l - lr
               ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"

           enddo
         enddo
         enddo

      else

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           lij =       inddm( ind_loop(1) , j             , k )
           lr  = lij - inddm( ind_loop(1) , jref - j      , k )
           lp  = lij - indven(ind_loop(1) , ind_loop(4)+1 , k )

!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
               li   = indbci(l-lij+ind_loop(1),  k )
               ldjr = l - lr
               ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank_SA.for"

           enddo
         enddo
         enddo
       endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

      jref = 2*ind_loop(3) - 1

      if(param_int(NEQ).eq.5) then

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           lij =       inddm( ind_loop(1) , j             , k )
           lr  = lij - inddm( ind_loop(1) , jref - j      , k )
           lp  = lij - indven(ind_loop(1) , ind_loop(3)   , k )
!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
             li   = indbci(l-lij+ind_loop(1),  k )
             ldjr = l - lr
             ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"

           enddo
         enddo
         enddo

      else

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           lij =       inddm( ind_loop(1) , j           , k )
           lr  = lij - inddm( ind_loop(1) , jref - j    , k )
           lp  = lij - indven(ind_loop(1) , ind_loop(3) , k )
!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
             li   = indbci(l-lij+ind_loop(1),  k )
             ldjr = l - lr
             ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank_SA.for"

           enddo
         enddo
         enddo

       endif !param_int(NEQ)


      ELSEIF (idir.eq.5) THEN

      kref = 2*ind_loop(6) + 1

      if(param_int(NEQ).eq.5) then

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           lij =       inddm( ind_loop(1) , j  , k             )
           lr  = lij - inddm( ind_loop(1) , j  , kref - k      )
           lp  = lij - indven(ind_loop(1) , j ,  ind_loop(6)+1 )
!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
             li   = indbci(l-lij+ind_loop(1),  j )
             ldjr = l - lr
             ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"

           enddo
         enddo
         enddo
      else
        do k = ind_loop(5), ind_loop(6)
        do j = ind_loop(3), ind_loop(4)
          lij =       inddm( ind_loop(1) , j  , k             )
          lr  = lij - inddm( ind_loop(1) , j  , kref - k      )
          lp  = lij - indven(ind_loop(1) , j ,  ind_loop(6)+1 )
       !DEC$ IVDEP
          do l = lij, lij + ind_loop(2) - ind_loop(1)
            li   = indbci(l-lij+ind_loop(1),  j )
            ldjr = l - lr
            ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank_SA.for"

          enddo
        enddo
        enddo

       endif !param_int(NEQ)

      ELSE

      kref = 2*ind_loop(5) - 1

      if(param_int(NEQ).eq.5) then


         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           lij =       inddm( ind_loop(1) , j  , k           )
           lr  = lij - inddm( ind_loop(1) , j  , kref - k    )
           lp  = lij - indven(ind_loop(1) , j ,  ind_loop(5) )
!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
             li   = indbci(l-lij+ind_loop(1),  j )
             ldjr = l - lr
             ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"
           enddo
         enddo
         enddo

      else

         do k = ind_loop(5), ind_loop(6)
         do j = ind_loop(3), ind_loop(4)
           lij =       inddm( ind_loop(1) , j  , k           )
           lr  = lij - inddm( ind_loop(1) , j  , kref - k    )
           lp  = lij - indven(ind_loop(1) , j ,  ind_loop(5) )
!DEC$ IVDEP
           do l = lij, lij + ind_loop(2) - ind_loop(1)
             li   = indbci(l-lij+ind_loop(1),  j )
             ldjr = l - lr
             ldp  = l - lp

#include     "FastS/BC/BCWallViscousIsothermal_firstrank_SA.for"

           enddo
         enddo
         enddo

       endif !param_int(NEQ)

      ENDIF !idir

      END
