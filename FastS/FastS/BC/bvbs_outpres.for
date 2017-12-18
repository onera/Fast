c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_outpres(idir,lrhs, neq_mtr, param_int, 
     &                        ind_loop,
     &                        param_real, c4,c5,c6,
     &                        ventijk, tijk, rop, pext,
     &                        size_data, inc_bc)
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
      INTEGER_E size_data, inc_bc 
      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E c4,c5,c6, param_real(0:*)

C...  pext contient les donnees CGNS pour la sortie subsonique
      REAL_E pext(size_data)

C Var local
      INTEGER_E l,lij,lr,lp,i,j,k,l1,l2,ic,jc,kc,kc_vent,ldp,lmtr,l0,
     & incj, inck, li, nd

      !INTEGER_E k_piv,nd

      REAL_E ro,u,v,w,t,nut,c6inv,c0,c1,c2,c3,roe,rue,rve,rwe,ete,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & qref1,qref2,qref3,qref4,qref5,r,p,c,snorm, nue, gam1,
     & gam1_1,gam,cvinv,gamm1_1,gamm1,
     & sni,tn,qn,qen,eigv1,eigv2,eigv3,eigv4,
     & eigv5,qvar1,qvar2,qvar3,qvar4,qvar5,svar1,svar2,svar3,
     & svar4,svar5,rvar1,rvar2,rvar3,rvar4,rvar5,
     & roext,ruext,rvext,rwext,etext,roint,ruint,rvint,rwint,etint,
     & s_1,roi,rui,rvi,rwi,eti, tnx,tny, tnz, ri, roinv, sn, roe_inv
      REAL_E pi,ci,roci,dqn,rog,ug,vg,wg,Tg,pg,rgp,roc0 
      REAL_E ui,vi,wi,Ti,ro0,u0,v0,w0,T0,p0

C...  Tableau de travail local
      !REAL_E pext(size_data)
      
      INTEGER_E indbci
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"
      indbci(j_1,k_1) = 1 + (j_1-1) + (k_1-1)*inc_bc

      !do  nd = 1, size_data
      !   pext(nd) = data_pres(nd)
      !enddo
      !k_piv = data_pres(size_data+1)
         
c......determine la forme des tableuz metrique en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      snorm =-1.
      if(mod(idir,2).eq.0) snorm = 1.
      
      gam     = param_real(GAMMA)
      gam1    = param_real(GAMMA) - 1.
      gam1_1  = 1. / gam1
      rgp     = param_real(CVINF)*gam1
      
      ! semble inutile car pas d extrapolation ordre 3 dans la 2eme couche de ghost cells 
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
               li   = indbci(j,  k )
#include       "FastS/BC/BCOutpres_firstrank.for"

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
               li   = indbci(j,  k )
#include       "FastS/BC/BCOutpres_firstrank_SA.for"

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
               li   = indbci(j,  k )
#include       "FastS/BC/BCOutpres_firstrank.for"

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
               li   = indbci(j,  k )
#include       "FastS/BC/BCOutpres_firstrank_SA.for"
 

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
                  li   = indbci(i,  k )
#include          "FastS/BC/BCOutpres_firstrank.for"
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
                  li   = indbci(i,  k )
#include          "FastS/BC/BCOutpres_firstrank_SA.for"
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
                  li   = indbci(i,  k )
#include          "FastS/BC/BCOutpres_firstrank.for"
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
                  li   = indbci(i,  k )
#include          "FastS/BC/BCOutpres_firstrank_SA.for"
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
                  li   = indbci(i,  j )
#include        "FastS/BC/BCOutpres_firstrank.for"
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
                  li   = indbci(i,  j )
#include        "FastS/BC/BCOutpres_firstrank_SA.for"
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
#include        "FastS/BC/BCOutpres_firstrank_SA.for"
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
                  li   = indbci(i,  j )
#include        "FastS/BC/BCOutpres_firstrank.for"
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
                  li   = indbci(i,  j )
#include        "FastS/BC/BCOutpres_firstrank_SA.for"
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
