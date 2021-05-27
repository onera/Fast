c***********************************************************************
c    $Elsa-like injection massflow rate BC
c    $Author: RomainParis$
c***********************************************************************
      subroutine bvbs_injmfr(idir,lrhs, neq_mtr, param_int,
     &                        ind_loop,
     &                        param_real, c4,c5,c6,
     &                        ventijk, tijk, rop, d0x,
     &                        d0y, d0z, qp, ha, nue,
     &                        size_data, inc_bc)
c***********************************************************************
c put info here
c***********************************************************************
      implicit none


#include "FastS/param_solver.h"

      INTEGER_E idir,lrhs, neq_mtr, ind_loop(6), param_int(0:*)
      INTEGER_E size_data, inc_bc(3)
      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E c4,c5,c6, param_real(0:*)

C...  Valeurs des grandeurs thermo imposees
      REAL_E qp(size_data)                      ! Mass flow per area unit (B.C) defined on the interfaces
      REAL_E ha(size_data)                      ! Total enthalpy
      REAL_E d0x(size_data)                     ! x-component of the flow direction (B.C)
      REAL_E d0y(size_data)                     ! y-component of the flow direction (B.C)
      REAL_E d0z(size_data)                     ! z-component of the flow direction (B.C)
      REAL_E nue(size_data)                     ! turbulent viscosity

C Var local
      INTEGER_E l,lij,lr,lp,i,j,k,l1,l2,ic,jc,kc,kc_vent,ldp,lmtr,l0,
     & incj, inck, li, nd

      !INTEGER_E k_piv,nd

      REAL_E ro,u,v,w,t,nut,c6inv,c0,c1,c2,c3,roe,rue,rve,rwe,ete,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & qref1,qref2,qref3,qref4,qref5,r,p,c,snorm, gam1,
     & gam1_1,gam,cvinv,gamm1_1,gamm1,gam_gam1,
     & sni,tn,qen,eigv1,eigv2,eigv3,eigv4,
     & eigv5,qvar1,qvar2,qvar3,qvar4,qvar5,svar1,svar2,svar3,
     & svar4,svar5,rvar1,rvar2,rvar3,rvar4,rvar5,
     & roext,ruext,rvext,rwext,etext,roint,ruint,rvint,rwint,etint,
     & s_1,roi,rui,rvi,rwi,eti, tnx,tny, tnz, ri, roinv, sn, roe_inv
      REAL_E pi,ci,roci,dqn,rog,ug,vg,wg,Tg,pg,rgp,roc0
      REAL_E ui,vi,wi,Ti,ro0,u0,v0,w0,T0,p0
      REAL_E ventx,venty,ventz

C computing "qn"
      REAL_E    coefa, coefb, coefc, delta                              ! Trinomial quantities
      REAL_E    usdn, usdn2                                             ! 1/(d.n), d velocity direction and n
                                                                        ! inward normal and 1/(d.n)**2
      REAL_E    qn                                                      ! Normal velocity in the boundary
      REAL_E    wn,f,df,residug,nitnwt,newtonmax,tolnewton,pa,Minf,b,dwn

      INTEGER_E indbci
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

      indbci(j_1,k_1) = 1 + (j_1-inc_bc(2)) + (k_1-inc_bc(3))*inc_bc(1)



c......determine la forme des tableaux metriques en fonction de la nature du domaine
      !Seule la valeur de k_vent et ck_vent a un sens dans cet appel
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      snorm =-1.
      if(mod(idir,2).eq.0) snorm = 1.

      gam     = param_real(GAMMA)
      gam1    = param_real(GAMMA) - 1.
      gam1_1  = 1. / gam1
      gam_gam1= gam / gam1
      rgp     = param_real(CVINF)*gam1

      IF (idir.eq.1) THEN                               ! idir -> looping on imin
          if(param_int(NEQ).eq.5) then                  ! param_int(NEQ) -> DNS/ 5 equations
             do  k = ind_loop(5), ind_loop(6)
             do  j = ind_loop(3), ind_loop(4)
               l    = inddm( ind_loop(2)    , j,  k )
               lmtr = indmtr(ind_loop(2)+1  , j,  k )
               ldp  = indven(ind_loop(2)+1  , j,  k )
               l1   = l + 1
               li   = indbci(j,  k )
#include       "FastS/BC/BCInjMFR_firstrank.for"
               l0   = l
               l2   = l + 2
               do i = ind_loop(1), ind_loop(2)-1
                  l    = inddm( i , j,  k )
#include          "FastS/BC/BC_nextrank.for"
               enddo	! i
             enddo	! j
             enddo	! k
          else						! param_int(NEQ) -> SA
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
               l    = inddm( ind_loop(2)    , j,  k )
               lmtr = indmtr(ind_loop(2)+1  , j,  k )
               ldp  = indven(ind_loop(2)+1  , j,  k )
               l1   = l  + 1
               li   = indbci(j,  k )
#include       "FastS/BC/BCInjMFR_firstrank_SA.for"
               l0   = l
               l2   = l + 2
               do  i = ind_loop(1), ind_loop(2)-1
                 l    = inddm( i , j,  k )
#include         "FastS/BC/BC_nextrank_SA.for"
               enddo	! i
             enddo	! j
             enddo	! k
          endif						! param_int(NEQ)


      ELSEIF (idir.eq.2) THEN				! idir -> looping on imax
          if(param_int(NEQ).eq.5) then			! param_int(NEQ) -> DNS/ 5 equations
            do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
               l    = inddm( ind_loop(1)    , j,  k )
               lmtr = indmtr(ind_loop(1)    , j,  k )
               ldp  = indven(ind_loop(1)    , j,  k )
               l1   = l  - 1
               li   = indbci(j,  k )
#include       "FastS/BC/BCInjMFR_firstrank.for"
               l0   = l
               l2   = l - 2
                do i = ind_loop(1)+1, ind_loop(2)
                 l    = inddm(  i , j, k )
#include         "FastS/BC/BC_nextrank.for"
                enddo  	! i
             enddo	! j
             enddo	! k
          else						! param_int(NEQ) -> SA
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
               l    = inddm( ind_loop(1)    , j,  k )
               lmtr = indmtr(ind_loop(1)    , j,  k )
               ldp  = indven(ind_loop(1)    , j,  k )
               l1   = l  - 1
               li   = indbci(j,  k )
#include       "FastS/BC/BCInjMFR_firstrank_SA.for"
               l0   = l
               l2   = l - 2
                do i = ind_loop(1)+1, ind_loop(2)
                 l    = inddm(  i , j, k )
#include        "FastS/BC/BC_nextrank_SA.for"
                enddo  	! i
             enddo	! j
             enddo	! k
          endif						! param_int(NEQ)


      ELSEIF (idir.eq.3) THEN				! idir -> looping on jmin
          incj = param_int(NIJK)			! increment on j
          if(param_int(NEQ).eq.5) then			! param_int(NEQ) -> DNS/ 5 equations
             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2)
                  l    = inddm( i ,  ind_loop(4)    , k )
                  lmtr = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp  = indven(i ,  ind_loop(4) +1 , k )
                  l1   = l +   incj
                  li   = indbci(i,  k )
#include          "FastS/BC/BCInjMFR_firstrank.for"
                enddo ! i
                do  j = ind_loop(3), ind_loop(4)-1
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2)
                     l    = inddm( i ,    j            , k )
                     l0   = inddm( i ,  ind_loop(4)    , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_nextrank.for"
                  enddo	! i
                enddo	! j
             enddo 	! k
          else						! param_int(NEQ) -> SA
             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2)
                  l    = inddm( i ,  ind_loop(4)    , k )
                  lmtr = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp  = indven(i ,  ind_loop(4) +1 , k )
                  l1   = l +   incj
                  li   = indbci(i,  k )
#include          "FastS/BC/BCInjMFR_firstrank_SA.for"
                enddo ! i
                do  j = ind_loop(3), ind_loop(4)-1
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2)
                     l    = inddm( i ,    j            , k )
                     l0   = inddm( i ,  ind_loop(4)    , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_nextrank_SA.for"
                  enddo	! i
                enddo 	! j
             enddo 	! k
          endif						! param_int(NEQ)


      ELSEIF (idir.eq.4) THEN				! idir -> looping on jmax
          incj = -param_int(NIJK)			! increment on j
          if(param_int(NEQ).eq.5) then			! param_int(NEQ) -> DNS/ 5 equations
             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2)
                  l    = inddm( i ,  ind_loop(3)  , k )
                  lmtr = indmtr(i ,  ind_loop(3)  , k )
                  ldp  = indven(i ,  ind_loop(3)  , k )
                  l1   = l +   incj
                  li   = indbci(i,  k )
#include          "FastS/BC/BCInjMFR_firstrank.for"
                enddo ! i
                do  j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2)
                     l    = inddm( i ,    j          , k )
                     l0   = inddm( i ,  ind_loop(3)  , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_nextrank.for"
                  enddo	! i
                enddo 	! j
             enddo 	! k
          else						! param_int(NEQ) -> SA
             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2)
                  l    = inddm( i ,  ind_loop(3)  , k )
                  lmtr = indmtr(i ,  ind_loop(3)  , k )
                  ldp  = indven(i ,  ind_loop(3)  , k )
                  l1   = l +   incj
                  li   = indbci(i,  k )
#include          "FastS/BC/BCInjMFR_firstrank_SA.for"
                enddo ! i
                do  j = ind_loop(3)+1, ind_loop(4)
!DEC$ IVDEP
                  do i = ind_loop(1), ind_loop(2)
                     l    = inddm( i ,    j          , k )
                     l0   = inddm( i ,  ind_loop(3)  , k )
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_nextrank_SA.for"
                  enddo	! i
                enddo 	! j
             enddo 	! k
          endif						! param_int(NEQ)


      ELSEIF (idir.eq.5) THEN				! idir -> looping on kmin
          inck = param_int(NIJK)*param_int(NIJK+1)	! increment on k
          if(param_int(NEQ).eq.5) then			! param_int(NEQ) -> DNS/ 5 equations
             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)
                  l    = inddm( i , j,  ind_loop(6)   )
                  lmtr = indmtr(i , j,  ind_loop(6)+1 )
                  ldp  = indven(i , j,  ind_loop(6)+1 )
                  l1   = l +   inck
                  li   = indbci(i,  j )
#include          "FastS/BC/BCInjMFR_firstrank.for"
               enddo ! i
             enddo   ! j
             do  k = ind_loop(5), ind_loop(6)-1
               do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
		   do  i = ind_loop(1), ind_loop(2)

		     l    = inddm( i , j,  k          )
		     l0   = inddm( i , j, ind_loop(6) )
		     l1   = l0 +   inck
		     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank.for"
		   enddo ! i
              enddo	 ! j
             enddo 	 ! k
          else						! param_int(NEQ) -> SA
             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)
                  l    = inddm( i , j,  ind_loop(6)   )
                  lmtr = indmtr(i , j,  ind_loop(6)+1 )
                  ldp  = indven(i , j,  ind_loop(6)+1 )
                  l1   = l +   inck
                  li   = indbci(i,  j )
#include        "FastS/BC/BCInjMFR_firstrank_SA.for"
               enddo ! i
             enddo   ! j
             do  k = ind_loop(5), ind_loop(6)-1
                do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
                   do  i = ind_loop(1), ind_loop(2)

                     l    = inddm( i , j,  k          )
                     l0   = inddm( i , j, ind_loop(6) )
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank_SA.for"
                   enddo ! i
               enddo	 ! j
             enddo	 ! k
          endif						! param_int(NEQ)


      ELSE						! idir -> looping on kmax
          inck = -param_int(NIJK)*param_int(NIJK+1)	! increment on k
          if(param_int(NEQ).eq.5) then			! param_int(NEQ) -> DNS/ 5 equations
             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)
                  l    = inddm( i , j,  ind_loop(5)   )
                  lmtr = indmtr(i , j,  ind_loop(5)   )
                  ldp  = indven(i , j,  ind_loop(5)   )
                  l1   = l +   inck
                  li   = indbci(i,  j )
#include        "FastS/BC/BCInjMFR_firstrank.for"
               enddo ! i
             enddo   ! j
             do  k = ind_loop(5)+1, ind_loop(6)
                do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
                   do  i = ind_loop(1), ind_loop(2)
                     l    = inddm( i , j,  k          )
                     l0   = inddm( i , j, ind_loop(5) )
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank.for"
                   enddo ! i
               enddo	 ! j
             enddo	 ! k
          else						! param_int(NEQ) -> SA
             do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
               do  i = ind_loop(1), ind_loop(2)
                  l    = inddm( i , j,  ind_loop(5)   )
                  lmtr = indmtr(i , j,  ind_loop(5)   )
                  ldp  = indven(i , j,  ind_loop(5)   )
                  l1   = l +   inck
                  li   = indbci(i,  j )
#include        "FastS/BC/BCInjMFR_firstrank_SA.for"
               enddo ! i
             enddo   ! j
             do  k = ind_loop(5)+1, ind_loop(6)
                do  j = ind_loop(3), ind_loop(4)
!DEC$ IVDEP
                   do  i = ind_loop(1), ind_loop(2)
                     l    = inddm( i , j,  k          )
                     l0   = inddm( i , j, ind_loop(5) )
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank_SA.for"
                   enddo ! i
               enddo	 ! j
             enddo	 ! k
          endif						! param_int(NEQ)
      ENDIF						! idir
      END
