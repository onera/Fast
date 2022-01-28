c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine corr_debit_ibm(ndom,idir, neq_mtr,ithread, nthread_max,
     &                       param_int, size_fen, facelist,
     &                       rop, tijk, vol, celln, flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,idir, ithread, nthread_max, size_fen, neq_mtr
      INTEGER_E param_int(0:*)
      INTEGER_E facelist(*)

      REAL_E rop(param_int(NDIMDX), param_int(NEQ) )
      REAL_E celln(param_int(NDIMDX))
      REAL_E tijk( param_int(NDIMDX_MTR) , neq_mtr )
      REAL_E vol( param_int(NDIMDX_MTR) )
      REAL_E flux(7)

c  Var loc
      INTEGER_E f,l,l1,l2,lt,chunk, ideb,ifin,adr_flu,i,j,k,l0,
     & inc,flag,ic,jc,kc,kc_vent,inci,incj,inck
      REAL_E c4,c5,c6,c7,vcorr,vtg,qm1,qp1,vcorr0,vtg0,vmin,tcx,tcy,tcz,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,u0,ut,vt,wt,qn,surf,s_1,u,v,w,
     & sens

#include "FastS/formule_param.h"

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
      c7 = (c4-c5)/c6

      vmin = 1e-6

      vtg      = (-c4*flux(2) -  c6*flux(3))/c5
      vtg0     = (-c4*vtg     -  c5*flux(2))/c6

      vcorr    = (flux(1)-vtg )/flux(4)
      vcorr0   = (flux(7)-vtg0)/flux(4)

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

            !on calcule l'inverse de la surface de la facette
            surf = sqrt(tcx*tcx + tcy*tcy +tcz*tcz)
            surf = max(surf,1e-30)
            s_1  = 1./surf

            tcx = tcx*s_1
            tcy = tcy*s_1
            tcz = tcz*s_1

            !! Qn_L= 0 --> rop(l)  
            u = c5*rop(l, 2) + c4*rop(l1,2) + c6*rop(l2,2)
            v = c5*rop(l, 3) + c4*rop(l1,3) + c6*rop(l2,3)
            w = c5*rop(l, 4) + c4*rop(l1,4) + c6*rop(l2,4)

            qn = (u*tcx+v*tcy+w*tcz)
            !qn = u*tcx

            ut = u-qn*tcx
            vt = v-qn*tcy
            wt = w-qn*tcz

            rop(l , 2)=(ut -c4*rop(l1,2) - c6*rop(l2,2))/c5
            rop(l , 3)=(vt -c4*rop(l1,3) - c6*rop(l2,3))/c5
            rop(l , 4)=(wt -c4*rop(l1,4) - c6*rop(l2,4))/c5

            !! Qn_R= 0 --> rop(l0)  
            u = c4*rop(l, 2) + c5*rop(l1,2) + c6*rop(l0,2)
            v = c4*rop(l, 3) + c5*rop(l1,3) + c6*rop(l0,3)
            w = c4*rop(l, 4) + c5*rop(l1,4) + c6*rop(l0,4)

            qn = (u*tcx+v*tcy+w*tcz)
            !qn = u*tcx

            ut = u-qn*tcx
            vt = v-qn*tcy
            wt = w-qn*tcz

            rop(l0 , 2)=(ut -c5*rop(l1,2) - c4*rop(l,2))/c6
            rop(l0 , 3)=(vt -c5*rop(l1,3) - c4*rop(l,3))/c6
            rop(l0 , 4)=(wt -c5*rop(l1,4) - c4*rop(l,4))/c6

            u0      = rop(l0, 2)*tcx + rop(l0, 3)*tcy + rop(l0, 4)*tcz
         enddo

      ELSEIF (param_int(ITYPZONE).eq.2) THEN

         do  f = ideb, ifin

            l     = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc
            lt    = 1
            
            tcx  = tijk(lt,ic)*ci_mtr
            tcy  = tijk(lt,jc)*cj_mtr
            tcz  = tijk(lt,kc)*ck_mtr

            !on calcule l'inverse de la surface de la facette
            surf = sqrt(tcx*tcx + tcy*tcy +tcz*tcz)
            surf = max(surf,1e-30)
            s_1  = 1./surf

            tcx = tcx*s_1
            tcy = tcy*s_1
            tcz = tcz*s_1

            !! Qn_L= 0 --> rop(l)  
            u = c5*rop(l, 2) + c4*rop(l1,2) + c6*rop(l2,2)
            v = c5*rop(l, 3) + c4*rop(l1,3) + c6*rop(l2,3)
            w = c5*rop(l, 4) + c4*rop(l1,4) + c6*rop(l2,4)

            qn = (u*tcx+v*tcy+w*tcz)

            ut = u-qn*tcx
            vt = v-qn*tcy
            wt = w-qn*tcz

            rop(l , 2)=(ut -c4*rop(l1,2) - c6*rop(l2,2))/c5
            rop(l , 3)=(vt -c4*rop(l1,3) - c6*rop(l2,3))/c5
            rop(l , 4)=(wt -c4*rop(l1,4) - c6*rop(l2,4))/c5

            !! Qn_R= 0 --> rop(l0)  
            u = c4*rop(l, 2) + c5*rop(l1,2) + c6*rop(l0,2)
            v = c4*rop(l, 3) + c5*rop(l1,3) + c6*rop(l0,3)
            w = c4*rop(l, 4) + c5*rop(l1,4) + c6*rop(l0,4)

            qn = (u*tcx+v*tcy+w*tcz)

            ut = u-qn*tcx
            vt = v-qn*tcy
            wt = w-qn*tcz

            rop(l0 , 2)=(ut -c5*rop(l1,2) - c4*rop(l,2))/c6
            rop(l0 , 3)=(vt -c5*rop(l1,3) - c4*rop(l,3))/c6
            rop(l0 , 4)=(wt -c5*rop(l1,4) - c4*rop(l,4))/c6

            u0      = rop(l0, 2)*tcx + rop(l0, 3)*tcy + rop(l0, 4)*tcz
         enddo

      ELSEIF (param_int(ITYPZONE).eq.3) THEN

         do  f = ideb, ifin

            l     = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc
            lt    = 1
            
            tcx  = tijk(lt,ic)*ci_mtr
            tcy  = tijk(lt,jc)*cj_mtr

            !on calcule l'inverse de la surface de la facette
            surf = sqrt(tcx*tcx + tcy*tcy)
            surf = max(surf,1e-30)
            s_1  = 1./surf

            tcx = tcx*s_1
            tcy = tcy*s_1

            !! Qn_L= 0 --> rop(l)  
            u = c5*rop(l, 2) + c4*rop(l1,2) + c6*rop(l2,2)
            v = c5*rop(l, 3) + c4*rop(l1,3) + c6*rop(l2,3)

            qn = (u*tcx+v*tcy)

            ut = u-qn*tcx
            vt = v-qn*tcy

            rop(l , 2)=(ut -c4*rop(l1,2) - c6*rop(l2,2))/c5
            rop(l , 3)=(vt -c4*rop(l1,3) - c6*rop(l2,3))/c5

            !! Qn_R= 0 --> rop(l0)  
            u = c4*rop(l, 2) + c5*rop(l1,2) + c6*rop(l0,2)
            v = c4*rop(l, 3) + c5*rop(l1,3) + c6*rop(l0,3)

            qn = (u*tcx+v*tcy)

            ut = u-qn*tcx
            vt = v-qn*tcy

            rop(l0 , 2)=(ut -c5*rop(l1,2) - c4*rop(l,2))/c6
            rop(l0 , 3)=(vt -c5*rop(l1,3) - c4*rop(l,3))/c6

            u0      = rop(l0, 2)*tcx + rop(l0, 3)*tcy
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
