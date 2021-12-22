c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_debit_ibm(ndom,idir, neq_mtr,
     &                       ithread, nthread_max,nitcfg,
     &                       param_int, param_real, size_fen,facelist,
     &                       rop, tijk, vol, flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,idir, ithread, nthread_max,nitcfg,size_fen,neq_mtr
      INTEGER_E param_int(0:*)
      INTEGER_E facelist(*)

      REAL_E param_real(0:*)
      REAL_E rop(param_int(NDIMDX), param_int(NEQ) )
      REAL_E tijk( param_int(NDIMDX_MTR) , neq_mtr )

      REAL_E vol( param_int(NDIMDX_MTR) )
      REAL_E flux(7)

c  Var loc
      INTEGER_E f,l,l0,l1,l2,lt,chunk, ideb,ifin,inc,
     & i,j,k,flag,ic,jc,kc,kc_vent,inci,incj,inck

      REAL_E c4,c5,c6,vint,rgp,sens,u,u1,u2,u0,tcx,tcy,tcz,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,qn,v,w

#include "FastS/formule_param.h"

      !if(idir.ne.3) return

      !write(*,*)'fen',ind_fen
      !write(*,*)'lop',ind_loop
      !write(*,*)'inbc',ndom,idir,facelist(1)

      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)



      chunk = size_fen/nthread_max
      ideb  = (ithread-1)*chunk+1
      ifin  = ithread*chunk 
      if(ithread.eq.nthread_max) ifin = size_fen

c      if (ithread.eq.1.and.ndom.ge.162)
c     &  write(*,*)'inbc',ndom,idir,facelist(1),size_fen,nitcfg
c      write(*,'(a,3i5,a,4i5)')
c     & 'ideb,ifin',ideb,ifin,ifin-ideb+1,'idir,fam,ndom',
c     & idir,nitcfg,ndom,facelist(1)

      if(idir.eq.1) then
        inc = 1
        sens= 1
        flag= 1
      elseif(idir.eq.2) then
        inc = -1
        sens= -1
        flag=  0
      elseif(idir.eq.3) then
        inc = param_int(NIJK)
        sens= 1
        flag= 1
      elseif(idir.eq.4) then
        inc = -param_int(NIJK)
        sens= -1
        flag=  0
      elseif(idir.eq.5) then
        inc = param_int(NIJK)*param_int(NIJK+1)
        sens= 1
        flag= 1
      else
        inc = -param_int(NIJK)*param_int(NIJK+1)
        sens= -1
        flag=  0
      endif

      inci =1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)


      c4 =  5. / 6.
      c5 =  2. / 6.
      c6 = -1. / 6.

      rgp= param_real(CVINF)*(param_real(GAMMA)-1.)  !Cv(gama-1)= R (gas parfait)

      IF (param_int(ITYPZONE).eq.0) THEN

         do  f = ideb, ifin

            l     = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc
            lt    = l + inc*flag
            
            tcx  = tijk(lt,ic)*ci_mtr
            tcy  = tijk(lt,jc)*cj_mtr
            tcz  = tijk(lt,kc)*ck_mtr

            u = rop(l , 2)*tcx + rop(l , 3)*tcy + rop(l , 4)*tcz
            u1= rop(l1, 2)*tcx + rop(l1, 3)*tcy + rop(l1, 4)*tcz
            u2= rop(l2, 2)*tcx + rop(l2, 3)*tcy + rop(l2, 4)*tcz
            u0= rop(l0, 2)*tcx + rop(l0, 3)*tcy + rop(l0, 4)*tcz

            flux(1) = flux(1) +  u*sens
            flux(2) = flux(2) + u1*sens
            flux(3) = flux(3) + u2*sens
            flux(7) = flux(7) + u0*sens
            flux(4) = flux(4) +  tcx + tcy + tcz

            !u = c5*rop(l, 2) + c4*rop(l1,2) + c6*rop(l2,2)
            !v = c5*rop(l, 3) + c4*rop(l1,3) + c6*rop(l2,3)
            !w = c5*rop(l, 4) + c4*rop(l1,4) + c6*rop(l2,4)

            !k = l/inck
            !j = (l-k*inck)/incj
            !i = l-k*inck - j*incj
            !qn = (u*tcx+v*tcy+w*tcz)

            !! Qn_R= 0 --> rop(l0)  
            !u = c4*rop(l, 2) + c5*rop(l1,2) + c6*rop(l0,2)
            !v = c4*rop(l, 3) + c5*rop(l1,3) + c6*rop(l0,3)
            !w = c4*rop(l, 4) + c5*rop(l1,4) + c6*rop(l0,4)
            !qn = (u*tcx+v*tcy+w*tcz)

         enddo

      elseif  (param_int(ITYPZONE).eq.2) THEN !3Dcart

         do  f = ideb, ifin

            l     = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc
            lt    = 1
            
            tcx  = tijk(lt,ic)*ci_mtr
            tcy  = tijk(lt,jc)*cj_mtr
            tcz  = tijk(lt,kc)*ck_mtr

            u = rop(l , 2)*tcx + rop(l , 3)*tcy + rop(l , 4)*tcz
            u1= rop(l1, 2)*tcx + rop(l1, 3)*tcy + rop(l1, 4)*tcz
            u2= rop(l2, 2)*tcx + rop(l2, 3)*tcy + rop(l2, 4)*tcz
            u0= rop(l0, 2)*tcx + rop(l0, 3)*tcy + rop(l0, 4)*tcz

            flux(1) = flux(1) +  u*sens
            flux(2) = flux(2) + u1*sens
            flux(3) = flux(3) + u2*sens
            flux(7) = flux(7) + u0*sens
            flux(4) = flux(4) +  tcx + tcy + tcz
         enddo

      elseif  (param_int(ITYPZONE).eq.3) THEN !2Dcart

         do  f = ideb, ifin

            l     = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc
            lt    = 1
            
            tcx  = tijk(lt,ic)*ci_mtr
            tcy  = tijk(lt,jc)*cj_mtr

            u = c5*rop(l, 2) + c4*rop(l1,2) + c6*rop(l2,2)
            v = c5*rop(l, 3) + c4*rop(l1,3) + c6*rop(l2,3)

            k = l/inck
            j = (l-k*inck)/incj
            i = l-k*inck - j*incj
            qn = (u*tcx+v*tcy)

            flux(1) = flux(1) + (rop(l , 2)*tcx +rop(l , 3)*tcy)*sens
            flux(2) = flux(2) + (rop(l1, 2)*tcx +rop(l1, 3)*tcy)*sens
            flux(3) = flux(3) + (rop(l2, 2)*tcx +rop(l2, 3)*tcy)*sens
            flux(7) = flux(7) + (rop(l0, 2)*tcx +rop(l0, 3)*tcy)*sens
            flux(4) = flux(4) +  tcx + tcy

            !flux(1) = flux(1) +  rop(l , 2)*tijk(lt,1)*sens
            !flux(2) = flux(2) +  rop(l1, 2)*tijk(lt,1)*sens
            !flux(3) = flux(3) +  rop(l2, 2)*tijk(lt,1)*sens
            !flux(4) = flux(4) +  tijk(lt,1)
            !flux(7) = flux(7) +  rop(l0, 2)*tijk(lt,1)*sens
         enddo

      else !2Dcart
       if (ithread.eq.1) then
         write(*,*)'cp_debit pas code en 3Dhomo'
         stop
       endif
      endif

c      IF (idir.eq.1) THEN

c         do  f = ideb, ifin
c
c            l     = facelist( f)
c            l1    = facelist( f) + inci
c            l2    = facelist( f) + inci*2
c            l0    = facelist( f) - inci
c            lt    = 1
c            
c            flux(1) = flux(1) +  rop(l , 2)*tijk(lt,1)
c            flux(2) = flux(2) +  rop(l1, 2)*tijk(lt,1)
cc            flux(3) = flux(3) +  rop(l2, 2)*tijk(lt,1)
c            flux(4) = flux(4) +  tijk(lt,1)
c            flux(7) = flux(7) +  rop(l0, 2)*tijk(lt,1)
c            !flux(5) = flux(5) +  ti(lt,1)*rop(l,1)*rop(l, 5)*rgp
c
c            !k = l/inck
c            !j = (l-k*inck)/incj
c            !i = l-k*inck - j*incj
c            !l1= inddm(i-1,j-1,1)
c c        enddo
cc
c      ELSEIF (idir.eq.2) THEN
c
c         do  f = ideb, ifin
c
c            l     = facelist( f)
c            l1    = facelist( f) - inci
c            l2    = facelist( f) - inci*2
c            l0    = facelist( f) + inci
c            lt    = 1
c
c            flux(1) = flux(1) -  rop(l , 2)*ti(lt,1)
c            flux(2) = flux(2) -  rop(l1, 2)*ti(lt,1)
c           flux(3) = flux(3) -  rop(l2, 2)*ti(lt,1)
c           flux(7) = flux(7) -  rop(l0 ,2)*ti(lt,1)
c           flux(4) = flux(4) +  ti(lt,1)
            !flux(5) = flux(5) -  ti(lt,1)*rop(l,1)*rop(l,5)*rgp
            
c      if(nitcfg.eq.1.and.ndom.eq.26) write(*,'(a,3f15.11,i5)')"t1M",
c     & rop(l0,2), rop(l0 ,2)*ti(lt,1), flux(7), l0
c     &  rop(l,2)*ti(lt,1), flux(1),flux(4),l,i,j,k,ndom

c         enddo
c            write(*,'(a,4f15.10)')'Vmoy',vmoy1/float(l),vmoy2/float(l),
c     &                           vmoy3/float(l),vmoy4/float(l)

c      ELSEIF (idir.eq.3) THEN

c        do  f = ideb, ifin
         !if (ithread.eq.1.and.ndom.eq.0)write(*,*)'flis',f, facelist(f)

c            l     = facelist( f)
c            l1    = facelist( f) + incj
c            l2    = facelist( f) + incj*2
c            l0    = facelist( f) - incj
c            lt    = 1
c
c            flux(1) = flux(1) +  rop(l , 3)*tj(lt,2)
c            flux(2) = flux(2) +  rop(l1, 3)*tj(lt,2)
c            flux(3) = flux(3) +  rop(l2, 3)*tj(lt,2)
c            flux(7) = flux(7) +  rop(l0, 3)*tj(lt,2)
c            flux(4) = flux(4) +  tj(lt,2)
c            flux(6) = flux(6) +  tj(lt,2)*rop(l,1)*rop(l,5)*rgp
c            !flux(6) = flux(6) +  tj(lt,2)
c
cc      if(nitcfg.eq.2) write(*,'(a,4f15.11,5i5)')"t3m",rop(l,3),
c     &  rop(l,3)*tj(lt,2), flux(1),flux(4),l,i,j,k,ndom
c      if(ndom.eq.0) write(*,'(a,3f15.11,5i5)')"t1",rop(l,3),
c     &  rop(l,3)*tj(lt,2), flux(1),l,i,j,k,ndom

c         enddo

c      ELSEIF (idir.eq.4) THEN
c
c         do  f = ideb, ifin
c
c            l     = facelist( f)
c            l1    = facelist( f) - incj
c            l2    = facelist( f) - incj*2
c            l0    = facelist( f) + incj
c            lt    = 1
c
c            flux(1) = flux(1) -  rop(l , 3)*tj(lt,2)
c            flux(2) = flux(2) -  rop(l1, 3)*tj(lt,2)
c            flux(3) = flux(3) -  rop(l2, 3)*tj(lt,2)
c            flux(7) = flux(7) -  rop(l0, 3)*tj(lt,2)
c            flux(4) = flux(4) +  tj(lt,2)
c            flux(6) = flux(6) -  tj(lt,2)*rop(l,1)*rop(l,5)*rgp

c      if(nitcfg.eq.2) write(*,'(a,4f15.11,5i5)')"t3m",rop(l,3),
c     &  rop(l,3)*tj(lt,2), flux(1),flux(4),l,i,j,k,ndom

c         enddo
 
c      ELSE

c       continue
c      ENDIF

c      vint   = (c5*flux(1) +c4*flux(2) + c6*flux(3))

      !if(ndom.eq.249.and.idir.eq.4.and.nitcfg.eq.3)
c      if(ndom.le.249.and.idir.le.4)
c     &  write(*,'(a,5f18.15,a,2i4)')
c     & 'flux paroi',flux(7), flux(1) ,flux(2) ,flux(3),flux(4)
c     & ,'ndom, idir=',ndom,idir

      end
