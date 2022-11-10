c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine corr_debit_ibm(ndom,idir, neq_mtr,ithread, nthread_max,
     &                       param_int, param_real,
     &                       ipass,fam, nitcfg, amor,size_fen,
     &                       facelist, rop, tijk, vol, celln, flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,idir, ithread, nthread_max, size_fen, neq_mtr
      INTEGER_E param_int(0:*)
      INTEGER_E facelist(*),fam, nitcfg, ipass

      REAL_E rop(param_int(NDIMDX), param_int(NEQ) )
      REAL_E celln(param_int(NDIMDX))
      REAL_E tijk( param_int(NDIMDX_MTR) , neq_mtr )
      REAL_E vol( param_int(NDIMDX_MTR) )
      REAL_E flux(7),amor(2), param_real(0:*)

c  Var loc
      INTEGER_E f,l,l1,l2,lt,chunk, ideb,ifin,adr_flu,i,j,k,l0,
     & inc,flagi,flagj,flagk,ic,jc,kc,kc_vent,inci,incj,inck,var_tg,
     & var_mtr
      REAL_E c4,c5,c6,c7,vcorr,vtg,qm1,qp1,vcorr0,vtg0,tcx,tcy,tcz,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,u0,ut,vt,wt,qn,surf,s_1,u,v,w,
     & sens,ec_old,ec0_old,cpinv,pold,pold0

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !if(idir.gt.3.and.ndom.eq.42) return
      !return
      !write(*,*)'fen',ind_fen
      !write(*,*)'lop',ind_loop

      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      chunk = size_fen/nthread_max
      ideb  = (ithread-1)*chunk+1
      ifin  = ithread*chunk 
      if(ithread.eq.nthread_max) ifin = size_fen
 
      !ideb = 1
      !ifin = size_fen

      if(idir.eq.1) then
        inc     = 1
        sens    = 1
        flagi   = 1
        flagj   = 0
        flagk   = 0
        var_tg  = 2
        var_mtr = ic
      elseif(idir.eq.2) then
        inc     = -1
        sens    = -1
        flagi   = 0
        flagj   = 0
        flagk   = 0
        var_tg  = 2
        var_mtr = ic
      elseif(idir.eq.3) then
        inc     = param_int(NIJK)
        sens= 1
        flagi   = 0
        flagj   = 1
        flagk   = 0
        var_tg  = 3
        var_mtr = jc
      elseif(idir.eq.4) then
        inc     = -param_int(NIJK)
        sens= -1
        flagi   = 0
        flagj   = 0
        flagk   = 0
        var_tg  = 3
        var_mtr = jc
      elseif(idir.eq.5) then
        inc     = param_int(NIJK)*param_int(NIJK+1)
        sens= 1
        flagi   = 0
        flagj   = 0
        flagk   = 1
        var_tg  = 4
        var_mtr = kc
      else
        inc     = -param_int(NIJK)*param_int(NIJK+1)
        sens= -1
        flagi   = 0
        flagj   = 0
        flagk   = 0
        var_tg  = 4
        var_mtr = kc
      endif

      inci =1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      c4 =  5. / 6.
      c5 =  2. / 6.
      c6 = -1. / 6.
      c7 = (c4-c5)/c6

      cpinv = 1./(param_real( GAMMA )*param_real( CVINF ))

      !vtg      = (-c4*flux(2)    -  c6*flux(3))/c5/amor(1)
      !vtg0     = (-c4*vtg*amor(1)-  c5*flux(2))/c6/amor(2)

      vtg      = (-c4*flux(2) -  c6*flux(3))/c5
      vtg0     = (-c4*vtg     -  c5*flux(2))/c6

      vcorr    = (flux(1)-vtg )/flux(4)
      vcorr0   = (flux(7)-vtg0)/flux(4)

      !if (nitcfg.gt.1) then
      !   vcorr = vcorr  + amor(1)
      !   vcorr0= vcorr0 + amor(2)
      !endif


      vcorr  = amor(1)
      vcorr0 = amor(2)




c      if(ithread.eq.1.and.fam.eq.0) write(*,*)"vcorr",vcorr,vcorr0,
c     &   ndom,nitcfg,idir

      IF (param_int(ITYPZONE).eq.0) THEN

         do  f = ideb, ifin

            l     = facelist( f)
            l0    = facelist( f) - inc


            k = (l-1)/inck -1
            j = (l-1-(k+1)*inck)/incj -1
            i = (l-1-(k+1)*inck-(j+1)*incj) -1

            lt    = indmtr(i+flagi,j+flagj,k+flagk)

            s_1 = 1./sqrt(tijk(lt,1)**2+tijk(lt,2)**2+tijk(lt,3)**2)
           
            !rop(l , var_tg)= rop(l , var_tg) - vcorr*sens
            !rop(l0, var_tg)= rop(l0, var_tg) - vcorr0*sens

            !write(*,*)"argh ", ndom, i,j,k,celln(l),celln(l0)
c      if(ithread.eq.1.and.fam.eq.0.and.ndom.eq.0.and.f.eq.ideb
c     & .and.idir.eq.2) then
c       write(*,'(a,5f15.8,i3)')"update", rop(l,2),rop(l0, 2),
c     &   vcorr*sens*tijk(lt,1)*s_1,rop(l,1),rop(l0,1),nitcfg
c      endif

            !rop(l,  1)= rop(l  , 1)+ vcorr
            !rop(l0, 1)= rop(l0 , 1)+ vcorr0
      
            ec_old  = (rop(l,2)**2+rop(l,3)**2+rop(l,4)**2)
            ec0_old = (rop(l0,2)**2+rop(l0,3)**2+rop(l0,4)**2)

            pold = rop(l,  1)*rop( l , 5)
            pold0= rop(l0,  1)*rop(l0, 5)

c            rop(l,  2)= rop(l , 2)+ vcorr*sens*tijk(lt,1)*s_1
c            rop(l,  3)= rop(l , 3)+ vcorr*sens*tijk(lt,2)*s_1
c            rop(l,  4)= rop(l , 4)+ vcorr*sens*tijk(lt,3)*s_1
c            rop(l0, 2)= rop(l0, 2)+ vcorr0*sens*tijk(lt,1)*s_1
c            rop(l0, 3)= rop(l0, 3)+ vcorr0*sens*tijk(lt,2)*s_1
c            rop(l0, 4)= rop(l0, 4)+ vcorr0*sens*tijk(lt,3)*s_1

c            rop(l,  5)= rop(l, 5)+ 0.5*cpinv*(ec_old- ( rop(l,2)**2
c     &                                                 +rop(l,3)**2
c     &                                                 +rop(l,4)**2))

c            rop(l0,  5)= rop(l0, 5)+ 0.5*cpinv*(ec0_old- ( rop(l0,2)**2
c     &                                                    +rop(l0,3)**2
c     &                                                   +rop(l0,4)**2))
 
           !rop(l,  1)  = pold/rop(l,  5)
           !rop(l0,  1) = pold0/rop(l0,  5)

           !rop(l,  3)= rop(l , 3)- vcorr*sens*tijk(lt,2)*s_1/rop(l,1)
            !rop(l,  4)= rop(l , 4)- vcorr*sens*tijk(lt,3)*s_1/rop(l,1)
            !rop(l0, 2)= rop(l0, 2)- vcorr0*sens*tijk(lt,1)*s_1/rop(l0,1)
            !rop(l0, 3)= rop(l0, 3)- vcorr0*sens*tijk(lt,2)*s_1/rop(l0,1)
            !rop(l0, 4)= rop(l0, 4)- vcorr0*sens*tijk(lt,3)*s_1/rop(l0,1)


            flux(5) = flux(5) -  rop(l0 , var_tg)*tijk(lt,var_mtr)
         enddo

      ELSEIF (param_int(ITYPZONE).ge.2) THEN !2dcart  et 3d cart)

         do  f = ideb, ifin

            l     = facelist( f)
            l0    = facelist( f) - inc
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2

            lt    = 1
            
            rop(l , var_tg)= rop(l , var_tg) - vcorr*sens/rop(l,1)
            rop(l0, var_tg)= rop(l0, var_tg) - vcorr0*sens/rop(l0,1)

            !rop(l , var_tg)= rop(l , var_tg) - vcorr*sens
            !rop(l0, var_tg)= rop(l0, var_tg) - vcorr0*sens
            !rop(l , var_tg)= -rop(l1 , var_tg) 
            !rop(l0, var_tg)= -rop(l2, var_tg) 
            !rop(l , 2)= -rop(l1 , 2) 
            !rop(l0, 2)= -rop(l2 , 2) 
            !rop(l , 3)= -rop(l1 , 3) 
            !rop(l0, 3)= -rop(l2 , 3) 
            !rop(l , 4)= -rop(l1 , 4) 
            !rop(l0, 4)= -rop(l2 , 4) 

            flux(5) = flux(5) -  rop(l0 , var_tg)*tijk(lt,var_mtr)
         enddo

      ELSE !2Dcart
       if (ithread.eq.1) then
         write(*,*)'corr_debit pas code en 3Dhomo'
         stop
       endif
      ENDIF

c      IF (idir.eq.1) THEN
c
c!CDIR novec
c         do  f = ideb, ifin
c
c            l     = facelist( f)
c            l1    = facelist( f) + inci
c            l2    = facelist( f) + inci*2
c            l3    = facelist( f) - inci
c            lt    = 1
c
c            !u_L=0
c            rop(l , 2)= (-c4*rop(l1,2) - c6*rop(l2,2))/c5
c            !u_R=-u_L
c            rop(l3, 2)=(-c6*rop(l2,2) -(c5+c4)*(rop(l1,2)+rop(l,2)))/c6
c            
c            !rop(l , 2)= (-c4*rop(l1,2) - c6*rop(l2,2))/c5
c            !rop(l3, 2)= (-c5*rop(l1,2) - c4*rop(l ,2))/c6
c
c         enddo
c
c      ELSEIF (idir.eq.2) THEN
c
c!CDIR novec
c         do  f = ideb, ifin
c
c            l     = facelist( f)                      !l
c            l1    = facelist( f) - inci               !nm
c            l2    = facelist( f) - inci*2             !nm2
c            l3    = facelist( f) + inci               !np
c            lt    = 1
c            !rop(l , 2)= (-c6-c4*0.5)*rop(l2,2)/(c5+c4*0.5)
c            !rop(l3, 2)= ((-c4-c5*0.5)*rop(l,2) - c5*0.5*rop(l2,2))/c6
c!
c            !rop(l1, 2)= 0.5*( rop(l2, 2) +  rop(l, 2) )
c            !rop(l1, 3)= 0.5*( rop(l2, 3) +  rop(l, 3) )
c
c            rop(l , 2)= (-c4*rop(l1,2) - c6*rop(l2,2))/c5
c            rop(l3, 2)=(-c6*rop(l2,2) -(c5+c4)*(rop(l1,2)+rop(l,2)))/c6
c
c            !rop(l , 2)= (-c4*rop(l1,2) - c6*rop(l2,2))/c5
c            !rop(l3, 2)= (-c5*rop(l1,2) - c4*rop(l ,2))/c6
c            !qm1=  c4*rop(l ,2)  +  c5*rop(nm , 2) +  c6*rop(np, 2)
c            !qm1=  c4*rop(l ,2)  +  c5*rop(l1 , 2) +  c6*rop(l3, 2)
c            !qp1=  c4*rop(nm, 2) +  c6*rop(nm2, 2) +  c5*rop(l , 2)
c            !qp1=  c4*rop(l1, 2) +  c6*rop(l2 , 2) +  c5*rop(l , 2)
c            !rop(l , 2)= rop(l ,2) + vcorr
c            !rop(l3, 2)= rop(l3,2) + vcorr0
c         enddo

c      ELSEIF (idir.eq.4) THEN
c
c         do  f = ideb, ifin

c            l     = facelist( f)
c            l1    = facelist( f) - incj
c            l2    = facelist( f) - incj*2
c            l3    = facelist( f) + incj
c            lt    = 1

c            rop(l , 3)=(-c4*rop(l1,3) - c6*rop(l2,3))/c5
c            rop(l3, 3)=(-c6*rop(l2,3) -(c5+c4)*(rop(l1,3)+rop(l,3)))/c6

c         enddo
 
c      ELSE

c       continue
c      ENDIF

      !if(ndom.eq.249.and.idir.eq.4)
c      if(ndom.le.249.and.idir.le.4)
c     &  write(*,'(a,f18.15,a,3f20.15,a,2i4)')
c     & 'verif surf',flux(5),'  flu(7),vtg, corr=',flux(7),vtg0,vcorr0
c     & ,'ndom, idir=',ndom,idir

      end
