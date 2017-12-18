c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_outflow(idir,lrhs, neq_mtr, param_int, ind_loop,
     &                         param_real, c4,c5,c6,
     &                         ventijk, tijk, rop, state)
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
      REAL_E state(param_int(NEQ))
      REAL_E c4,c5,c6, param_real(0:*)

C Var local
      INTEGER_E l,lij,lr,lp,i,j,k,l1,l2,ic,jc,kc,kc_vent,ldp,lmtr,l0,
     & incj, inck
      REAL_E ro,u,v,w,t,nut,c6inv,c0,c1,c2,c3,roe,rue,rve,rwe,ete,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & qref1,qref2,qref3,qref4,qref5,r,p,c,ci,snorm, nue, gamm1,
     & gamm1_1, cvinv,
     & sni,tn,qn,qen,eigv1,eigv2,eigv3,eigv4,
     & eigv5,qvar1,qvar2,qvar3,qvar4,qvar5,svar1,svar2,svar3,
     & svar4,svar5,rvar1,rvar2,rvar3,rvar4,rvar5,
     & roext,ruext,rvext,rwext,etext,roint,ruint,rvint,rwint,etint,
     $ s_1,roi,rui,rvi,rwi,eti, tnx,tny, ri, roinv, sn, roe_inv

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)


c      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop=', ind_loop

      snorm =-1.
      if(mod(idir,2).eq.0) snorm = 1.

      !!! premiere rangee = state
      !!! rangee suivante pour que : Phi_L + Phi_R = 2 * state

      roe=state(1)
      rue=state(2)
      rve=state(3)
      rwe=state(4)
      ete=state(5)

      if(param_int(NEQ).eq.6) nue = state(6)/roe

      roe_inv = 1./state(1)

      gamm1   = param_real(GAMMA) - 1.
      gamm1_1 = 1./gamm1
      cvinv   = 1./param_real(CVINF)

      c0   = 1./c6
      c1   =-(c4 + c5)*c0
      c2   =- c6*c0
      c3   = (2.- c5- c4)*c0

      IF (idir.eq.1) THEN

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
             do  j = ind_loop(3), ind_loop(4)

               l    = inddm( ind_loop(2)    , j,  k )
               lmtr = indmtr(ind_loop(2)+1  , j,  k )
               ldp  = indven(ind_loop(2)+1  , j,  k )
               l1   = l + 1

#include       "FastS/BC/BCOutflow_firstrank.for"

               l0   = l
               l2   = l + 2
               do i = ind_loop(1), ind_loop(2)-1

                  l    = inddm( i , j,  k ) 
#include          "FastS/BC/BC_extrap.for"
               enddo
             enddo
             enddo


          else

             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)

               l    = inddm( ind_loop(2)    , j,  k )
               lmtr = indmtr(ind_loop(2)+1  , j,  k )
               ldp  = indven(ind_loop(2)+1  , j,  k )
               l1   = l  + 1

#include       "FastS/BC/BCOutflow_firstrank_SA.for"

               l0   = l
               l2   = l + 2
               do  i = ind_loop(1), ind_loop(2)-1
 
                 l    = inddm( i , j,  k ) 
#include        "FastS/BC/BC_extrap_SA.for"
               enddo   
             enddo
             enddo
           endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

          if(param_int(NEQ).eq.5) then


             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)

               l    = inddm( ind_loop(1)    , j,  k )
               lmtr = indmtr(ind_loop(1)    , j,  k )
               ldp  = indven(ind_loop(1)    , j,  k )
               l1   = l  - 1
#include       "FastS/BC/BCOutflow_firstrank.for"

               l0   = l
               l2   = l - 2
                do i = ind_loop(1)+1, ind_loop(2)

                 l    = inddm(  i , j, k ) 
#include        "FastS/BC/BC_extrap.for"
                enddo   
             enddo
             enddo

          else

             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)

               l    = inddm( ind_loop(1)    , j,  k )
               lmtr = indmtr(ind_loop(1)    , j,  k )
               ldp  = indven(ind_loop(1)    , j,  k )
               l1   = l  - 1
#include       "FastS/BC/BCOutflow_firstrank_SA.for"
 

               l0   = l
               l2   = l - 2
                do i = ind_loop(1)+1, ind_loop(2)

                 l    = inddm(  i , j, k ) 
#include        "FastS/BC/BC_extrap_SA.for"
                enddo  
             enddo
             enddo
           endif !param_int(NEQ)


      ELSEIF (idir.eq.3) THEN

          incj = param_int(NIJK)

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l    = inddm( i ,  ind_loop(4)    , k )
                  lmtr = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp  = indven(i ,  ind_loop(4) +1 , k )
                  l1   = l +   incj
#include          "FastS/BC/BCOutflow_firstrank.for"
                enddo !i

                do  j = ind_loop(3), ind_loop(4)-1
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2) 

                     l    = inddm( i ,    j            , k )
                     l0   = inddm( i ,  ind_loop(4)    , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_extrap.for"
                  enddo!i
                enddo !j
             enddo !k

          else

             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l    = inddm( i ,  ind_loop(4)    , k )
                  lmtr = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp  = indven(i ,  ind_loop(4) +1 , k )
                  l1   = l +   incj
#include          "FastS/BC/BCOutflow_firstrank_SA.for"
                enddo !i

                do  j = ind_loop(3), ind_loop(4)-1
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2) 

                     l    = inddm( i ,    j            , k )
                     l0   = inddm( i ,  ind_loop(4)    , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_extrap_SA.for"
                  enddo!i
                enddo !j
             enddo !k

          endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

          incj = -param_int(NIJK)

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l    = inddm( i ,  ind_loop(3)  , k )
                  lmtr = indmtr(i ,  ind_loop(3)  , k )
                  ldp  = indven(i ,  ind_loop(3)  , k )
                  l1   = l +   incj
#include          "FastS/BC/BCOutflow_firstrank.for"
                enddo !i

                do  j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2) 

                     l    = inddm( i ,    j          , k )
                     l0   = inddm( i ,  ind_loop(3)  , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_extrap.for"
                  enddo!i
                enddo !j
             enddo !k

          else

             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                  l    = inddm( i ,  ind_loop(3)  , k )
                  lmtr = indmtr(i ,  ind_loop(3)  , k )
                  ldp  = indven(i ,  ind_loop(3)  , k )
                  l1   = l +   incj
#include          "FastS/BC/BCOutflow_firstrank_SA.for"
                enddo !i

                do  j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2) 

                     l    = inddm( i ,    j          , k )
                     l0   = inddm( i ,  ind_loop(3)  , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_extrap_SA.for"
                  enddo!i
                enddo !j
             enddo !k

          endif !param_int(NEQ)


      ELSEIF (idir.eq.5) THEN

          inck = param_int(NIJK)*param_int(NIJK+1)

          if(param_int(NEQ).eq.5) then

             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)

                  l    = inddm( i , j,  ind_loop(6)   )
                  lmtr = indmtr(i , j,  ind_loop(6)+1 )
                  ldp  = indven(i , j,  ind_loop(6)+1 )
                  l1   = l +   inck
#include        "FastS/BC/BCOutflow_firstrank.for"
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
#include             "FastS/BC/BC_extrap.for"
                   enddo
               enddo
             enddo

          else

             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)

                  l    = inddm( i , j,  ind_loop(6)   )
                  lmtr = indmtr(i , j,  ind_loop(6)+1 )
                  ldp  = indven(i , j,  ind_loop(6)+1 )
                  l1   = l +   inck
#include        "FastS/BC/BCOutflow_firstrank_SA.for"
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
#include             "FastS/BC/BC_extrap_SA.for"
                   enddo
               enddo
             enddo

          endif !param_int(NEQ)


      ELSE 

          inck = -param_int(NIJK)*param_int(NIJK+1)

          if(param_int(NEQ).eq.5) then

             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)

                  l    = inddm( i , j,  ind_loop(5)   )
                  lmtr = indmtr(i , j,  ind_loop(5)   )
                  ldp  = indven(i , j,  ind_loop(6)   )
                  l1   = l +   inck
#include        "FastS/BC/BCOutflow_firstrank.for"
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
#include             "FastS/BC/BC_extrap.for"
                   enddo
               enddo
             enddo

          else

             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)

                  l    = inddm( i , j,  ind_loop(5)   )
                  lmtr = indmtr(i , j,  ind_loop(5)   )
                  ldp  = indven(i , j,  ind_loop(6)   )
                  l1   = l +   inck
#include        "FastS/BC/BCOutflow_firstrank_SA.for"
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
#include             "FastS/BC/BC_extrap_SA.for"
                   enddo
               enddo
             enddo

          endif !param_int(NEQ)

      ENDIF !idir

      END
