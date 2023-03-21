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
     & incijk,ldjr2,lsample,lt,m,lmtr,iter,l0,ldjr3,ldjr4,ldnm,ldl,ldnp,
     & l1,l2,incj,inck,indbci,linput, lxyz, inck_x, incj_x

      REAL_E ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,tcx,tcy,tcz,
     & ventx,venty,ventz,r,u,v,w,ut,vt,wt,ua,va,wa,s_1,qn,dx,dy,dz,
     & nx, ny,nz,utau2,seuil,unorm1,qen,u_int,c0,c1,c2,c3,
     & dist,utau,unorm,surf,aa,yplus,g1,g2,g3,f,tp,fp,tauw,
     & Twall,uplus,pr,prtur,ratio,dt,pf,Tplus,qwall,cp,cond,vmoy_bc,
     & wall_ut,wall_vt,wall_wt,wall_int,co,si,rot1,rot2,ra,un,vn,wn,
     & umod,umoy,vmoy,wmoy,ufluc,vfluc,wfluc,ucut,vcut,wcut,rescale,
     & vmoy1, vmoy2, vmoy3, vmoy4,qm,qp,c7,vmoy2int
#include "FastS/formule_param.h"
#include "FastS/formule_xyz_param.h"
#include "FastS/formule_vent_param.h"
#include "FastS/formule_mtr_param.h"

      indbci(j_1,k_1) = 1 + (j_1-inc_bc(2)) + (k_1-inc_bc(3))*inc_bc(1)

      incj_x = param_int(NIJK_XYZ)
      inck_x = param_int(NIJK_XYZ)* param_int(NIJK_XYZ+1)
c      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
c      write(*,*)'nijk=', nijk
c      write(*,*)'loop vis=', ind_loop, inc_bc

      sample = param_int(WM_SAMPLING)

      ucut = 60.
      vcut = param_real(PSIROE)
      wcut = 30.

c      if (int( param_wmles(3) ) .eq.1) then
c       vmoy_bc = 0
c      else
       vmoy_bc = param_wmles(4)
c      endif

      !write(*,*)'vmoy_bc',vcut

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

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

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

          enddo
          enddo
       else
        ! a faire
       endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

          do k = ind_loop(5), ind_loop(6)
          do j = ind_loop(3), ind_loop(4)

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

          enddo
          enddo



      ELSEIF (idir.eq.3) THEN

          incj = param_int(NIJK)
          inck = param_int(NIJK)* param_int(NIJK+1)

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

                  lxyz = indcg( i ,  ind_loop(4) +4 , k )
                  lmtr = indmtr(i ,  ind_loop(4) +1 , k )
                  ldp  = indven(i ,  ind_loop(4) +1 , k )

                  linput = indbci(i,k) 

                  !ldjr   = l +   incj
                  !ldjr2  = l + 2*incj
                  lsample= ldjr + (sample-1)*incj
C#include          "FastS/BC/BCWallExchange.for"
#include          "FastS/BC/BCWallExchangeMoy.for"
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
                  enddo!i
                enddo !j
             enddo !k

            vmoy1 =0
            vmoy2 =0
            vmoy3 =0
            vmoy4 =0
             do  k = ind_loop(5), ind_loop(6)
!DEC$ IVDEP
                do i = ind_loop(1), ind_loop(2) 

                 ! l    = inddm( i ,  ind_loop(4)+2    , k )
                 !vmoy1 = vmoy1+rop(l,3)
                  l    = inddm( i ,  ind_loop(4)+3    , k )
                 vmoy2 = vmoy2+rop(l,3)
                  l    = inddm( i ,  ind_loop(4)+4    , k )
                 vmoy3 = vmoy3+rop(l,3)
                  l    = inddm( i ,  ind_loop(4)+5    , k )
                 vmoy4 = vmoy4+rop(l,3)
                enddo !i
                enddo !k

             l0 =(ind_loop(6)-ind_loop(5)+1)*(ind_loop(2)-ind_loop(1)+1)
c            !qm = (c4*vmoy3 + c5*vmoy2 + c6*vmoy4)/float(l0)

             !vmoy3 = 0.4*vmoy4

             !vmoy2int = (0.001*float(l0)-c4*vmoy3 -  c6*vmoy4)/c5
             vmoy2int = (0.001*float(l0)-c4*vmoy4*0.4 -  c6*vmoy4)/c5
c
             do  k = ind_loop(5), ind_loop(6)
                do i = ind_loop(1), ind_loop(2) 
c
c   cut3bis
c
c                 l       = inddm( i ,  ind_loop(4)+4    , k )
c                 linput  = indbci(i,k)
c                 rop(l,3)= rop(l,3)- vmoy_bc + 0.001
c   cut3ter
                 l       = inddm( i ,  ind_loop(4)+3    , k )
                 rop(l,3)= rop(l,3)- vmoy2/float(l0)+ vmoy2int/float(l0)

                 l       = inddm( i ,  ind_loop(4)+4    , k )
                 rop(l,3)= rop(l,3)- vmoy3/float(l0)+0.4*vmoy4/float(l0)
                enddo !i
                enddo !k

             do  k = ind_loop(5), ind_loop(6)
                do i = ind_loop(1), ind_loop(2),2 
c   cut3ter
                 l       = inddm( i ,  ind_loop(4)+3    , k )
                 rop(l,1)= 0.5*(rop(l+1,1)+ rop(l-1,1))

                 l       = inddm( i ,  ind_loop(4)+4    , k )
                 rop(l,1)= 0.5*(rop(l+1,1)+ rop(l-1,1))

                enddo !i
                enddo !k


             vmoy3 = 0.4*vmoy4

             c7 = (c4-c5)/c6
             do  k = ind_loop(5), ind_loop(6)
                do i = ind_loop(1), ind_loop(2) 
                  l     = inddm( i , ind_loop(4)+2    , k )
                  ldnm  = inddm( i , ind_loop(4)+3,  k )
                  ldl   = inddm( i , ind_loop(4)+4,  k )
                  ldnp  = inddm( i , ind_loop(4)+5,  k )
                 rop(l,1) = c7*( rop(ldl,1) - rop(ldnm,1) )+ rop(ldnp,1)
                 rop(l,2) = c7*( rop(ldl,2) - rop(ldnm,2) )+ rop(ldnp,2)
                 rop(l,3) = c7*( rop(ldl,3) - rop(ldnm,3) )+ rop(ldnp,3)
                 rop(l,4) = c7*( rop(ldl,4) - rop(ldnm,4) )+ rop(ldnp,4)
                 rop(l,5) = c7*( rop(ldl,5) - rop(ldnm,5) )+ rop(ldnp,5)
c                 !rop(l,1) = 2.*rop(ldnm,1) - rop(ldl,1)
c                 !rop(l,2) = 2.*rop(ldnm,2) - rop(ldl,2)
c                 !rop(l,3) = 2.*rop(ldnm,3) - rop(ldl,3)
c                 !rop(l,4) = 2.*rop(ldnm,4) - rop(ldl,4)
c                 !rop(l,5) = 2.*rop(ldnm,5) - rop(ldl,5)
                 vmoy1 = vmoy1+rop(l,3)
                enddo !i
                enddo !k

        !qm1=  c4*rop(l  +v1)  +  c5*rop(nm  +v1) +  c6*rop(np +v1)
        !qp1=  c4*rop(nm +v1)  +  c6*rop(nm2 +v1) +  c5*rop(l  +v1)

          l =(ind_loop(6)-ind_loop(5)+1)*(ind_loop(2)-ind_loop(1)+1)

          qm = (c4*vmoy3    + c5*vmoy2int + c6*vmoy4)/float(l)
          qp = (c4*vmoy2int + c5*vmoy3    + c6*vmoy1)/float(l)
c          !qm = (c4*vmoy3 + 0.*vmoy2 + c6*vmoy4)/float(l)
c          !qp = (0.*vmoy2 + c5*vmoy3 + 0.*vmoy1)/float(l)

c       if (ind_loop(5).eq.1.and.ind_loop(2).ge.780) 
c     & write(*,'(a,6f15.9)')'bcvmoy',vmoy1/int(l),vmoy2/float(l),
c     & vmoy3/float(l), vmoy4/float(l),qm,qp

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
