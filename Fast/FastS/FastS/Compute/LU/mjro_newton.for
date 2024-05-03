c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mjro_newton(ndom, nitcfg, neq,ndimdx, nijk, ijkv,
     &                       lSA, imtr, cv, visco, ratiom,ind_loop,
     &                       rop_1,rop,drodm)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Mise a jour de la sous-iteration P => P+1
c
c     VAL
c_V    3ADI,LU sous itere au 1er/2eme ordre en tps
c_V    Processeur domaine
c_V    Steady/Unsteady (Magnus)
c
c     I/O
c      ro
c***********************************************************************
      implicit  none


      INTEGER_E ndom,neq,ndimdx, lSA, imtr, nijk(5),ijkv(3),ind_loop(6),
     & nitcfg
      REAL_E rop_1(ndimdx,neq),rop(ndimdx,neq),drodm(ndimdx,neq)
      REAL_E visco(5),cv,ratiom

C var loc
      INTEGER_E incmax,i,j,k,l,lij,imin,jmin,kmin,imax,jmax,kmax
      REAL_E ro_old,u_old,v_old,w_old,t_old,roe_old,nu_old,anulam,
     & r_1,cvinv,cvinv2,rho,rhou,rhov,rhoe,temp01,cmus1,coesut

#include "FastS/formule.h"

      cvinv  = 1./cv   
      cvinv2 = 0.5*cvinv

      imin = max(1,ind_loop(1))
      jmin = max(1,ind_loop(3))
      kmin = max(1,ind_loop(5))
      imax = min(ijkv(1),ind_loop(2))
      jmax = min(ijkv(2),ind_loop(4))
      kmax = min(ijkv(3),ind_loop(6))

      If(lSA.eq.1) Then

      cmus1  = visco(5)
      temp01 = 1./visco(4)
      coesut = visco(3) * (1.+cmus1*temp01)

       if(imtr.eq.3) then

#ifndef E_SCALAR_COMPUTER
      incmax=nijk(1)*nijk(2)*nijk(5)
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
      do l=incmax+1,ndimdx-incmax
#else
      do j = jmin, jmax
        lij  = inddm(imin, j, 1)
        do l = lij, lij + imax - imin
CC attention rop et rop_1 identique sauf si nitcfg =1
#endif
        ro_old    = rop(l,1)
        u_old     = rop(l,2)
        v_old     = rop(l,3)
        t_old     = rop(l,5)
        nu_old    = rop(l,6)

        rop_1(l,1) = rop(l,1) + drodm(l,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm(l,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm(l,3))*r_1
        rop_1(l,4) = 0.

        roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                           +v_old*v_old))

        rop_1(l,5)     = ( roe_old + drodm(l,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3))


        anulam = coesut*sqrt(rop_1(l,5)*temp01)/(1.+cmus1/rop_1(l,5))
        anulam = anulam*r_1*ratiom 

        t_old  = (ro_old*nu_old + drodm(l,6))*r_1

        rop_1(l,6) = t_old
#ifndef E_SCALAR_COMPUTER
        enddo
#else
        enddo
        enddo
#endif

       else
#ifndef E_SCALAR_COMPUTER
      incmax=nijk(1)*nijk(2)*nijk(5)
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
      do l=incmax+1,ndimdx-incmax
#else
      do k = kmin, kmax
      do j = jmin, jmax
        lij  = inddm(imin , j, k)
        do l = lij, lij +  imax - imin
CC attention rop et rop_1 identique sauf si nitcfg =1
#endif
        ro_old     = rop(l,1)
        u_old      = rop(l,2)
        v_old      = rop(l,3)
        w_old      = rop(l,4)
        t_old      = rop(l,5)
        nu_old     = rop(l,6)

        rop_1(l,1) = rop(l,1) + drodm(l,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm(l,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm(l,3))*r_1
        rop_1(l,4) = (ro_old*w_old + drodm(l,4))*r_1

        roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                           +v_old*v_old
     &                                           +w_old*w_old))

        rop_1(l,5)     = ( roe_old +  drodm(l,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3)
     &                            +rop_1(l,4)*rop_1(l,4))

        anulam = coesut*sqrt(rop_1(l,5)*temp01)/(1.+cmus1/rop_1(l,5))
        anulam = anulam*r_1*ratiom 

        t_old  = (ro_old*nu_old + drodm(l,6))*r_1

        rop_1(l,6) = t_old

#ifndef E_SCALAR_COMPUTER
        enddo
#else
        enddo
        enddo
        enddo
#endif
       endif!2d/3d

      ELSE

       if(imtr.eq.3) then


#ifndef E_SCALAR_COMPUTER
      incmax=nijk(1)*nijk(2)*nijk(5)
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
      do l=incmax+1,ndimdx-incmax
#else
      do j = jmin, jmax
        lij  = inddm(imin, j, 1)
CDIR$ NOVEC
        do l = lij, lij + imax - imin
CC attention rop et rop_1 identique sauf si nitcfg =1
#endif
        ro_old    = rop(l,1)
        u_old     = rop(l,2)
        v_old     = rop(l,3)
        t_old     = rop(l,5)

        rop_1(l,1) = rop(l,1) + drodm(l,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm(l,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm(l,3))*r_1
        rop_1(l,4) = 0.

        roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                           +v_old*v_old))

        rop_1(l,5)     = ( roe_old + drodm(l,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3))
#ifndef E_SCALAR_COMPUTER
      enddo
#else
      enddo
      enddo
#endif

       else

#ifndef E_SCALAR_COMPUTER
      incmax=nijk(1)*nijk(2)*nijk(5)
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
      do l=incmax+1,ndimdx-incmax
#else
      do k = kmin, kmax
      do j = jmin, jmax
        lij  = inddm(imin , j, k)
        do l = lij, lij +  imax - imin
CC attention rop et rop_1 identique sauf si nitcfg =1
#endif
        ro_old     = rop(l,1)
        u_old      = rop(l,2)
        v_old      = rop(l,3)
        w_old      = rop(l,4)
        t_old      = rop(l,5)

        rop_1(l,1) = rop(l,1) + drodm(l,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm(l,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm(l,3))*r_1
        rop_1(l,4) = (ro_old*w_old + drodm(l,4))*r_1

        roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                           +v_old*v_old
     &                                           +w_old*w_old))

        rop_1(l,5)     = ( roe_old +  drodm(l,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3)
     &                            +rop_1(l,4)*rop_1(l,4))
#ifndef E_SCALAR_COMPUTER
       enddo
#else
        enddo
        enddo
        enddo
#endif
       endif!2d/3d
       ENDIF!neq=5 ou 6 

      end
