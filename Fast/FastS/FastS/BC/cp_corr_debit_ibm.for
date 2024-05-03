c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cp_corr_debit_ibm(ndom,idir, neq_mtr,
     &                       ithread, nthread_max,nitcfg,
     &                       param_int, param_real, size_fen,facelist,
     &                       rop, tijk, coe, drodm, flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,idir, ithread, nthread_max,nitcfg,size_fen,neq_mtr
      INTEGER_E param_int(0:*)
      INTEGER_E facelist(*)

      REAL_E param_real(0:*)
      REAL_E rop(param_int(NDIMDX)* param_int(NEQ) )
      REAL_E tijk( param_int(NDIMDX_MTR) , neq_mtr )

      REAL_E drodm(param_int(NDIMDX), param_int(NEQ) )
      REAL_E coe( param_int(NDIMDX), param_int(NEQ_COE) )
      REAL_E flux(7)

c  Var loc
      INTEGER_E f,l,l0,l1,l2,lt,chunk, ideb,ifin,inc,i,j,k,flag,ic,jc,
     & kc,kc_vent,inci,incj,inck,flagi,flagj,flagk,var_mtr,
     & v1,v2,v3,v4,v5,lc,nm,nm2,np,vslp,l4,deep,lb

      REAL_E c4,c5,c6,vint,rgp,sens,u1,u2,u0,tcx,tcy,tcz,
     & ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,qn,v,w,surf,rp,rm,
     & c1,c2,c3,qm,qp,
     & r1,h1,rou1,rov1,row1,r2,h2,rou2,rov2,row2,
     & gam,gam1,gam2,gam3,gam4,qn1,qn2,u,tdu,p1p2,roref,uref,tam,tam1,
     & qm1,qm2,qm3,qm4,qm5,qp1,qp2,qp3,qp4,qp5,
     & flu1,flu2,flu3,flu4,flu5,p1,p2,son,c,cp,
     & du,si,dv,dw,dp,dqn,s_1,nx,ny,nz,f1,f2,f3,f4,psiroe,r,q,h,r_1,
     & tmin_1,roff,delm,delp,delq,slp,slq,avmin,
     & flu1_old,flu2_old,flu3_old,flu4_old,flu5_old

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      avmin(c,r)=sign(1.,c)*max(0.,min(abs(c),sign(1.,c)*r))

      !-----Variables physiques
      gam    = param_real( GAMMA )
      rgp    = param_real( CVINF )*(gam-1.)  !Cv(gama-1)= R (gas parfait)
      gam1   = gam/(gam-1.)
      gam2   = 1./gam
      gam3    = gam1/ param_real( PRANDT )*rgp
      gam4    = gam1/ param_real( PRANDT_TUR )*rgp

      cp      =  param_real( CVINF )* param_real( GAMMA )
      roref= param_real( ROINF)
      uref = param_real( VINF )

      psiroe= param_real( PSIROE )
      tmin_1= 100./param_real( TINF )!!si T< 0.01Tinf, alors limiteur null

      c1     = 0.02*uref         ! modif suite chant metrique et suppression tc dans flux final
      c2     = 0.02/(uref*roref) ! modif suite chant metrique et suppression tc dans flux final
      c3     = -2.

      !    roff MUSCL
      c6     = 1./6.
      c4     = 5.*c6
      c5     = 2.*c6
      c6     =-1.*c6


      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

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
        inc     = 1
        sens    = 1
        flagi   = 1
        flagj   = 0
        flagk   = 0
        var_mtr = ic
      elseif(idir.eq.2) then
        inc     = -1
        sens    = -1
        flagi   = 0
        flagj   = 0
        flagk   = 0
        var_mtr = ic
      elseif(idir.eq.3) then
        inc     = param_int(NIJK)
        sens    = 1
        flagi   = 0
        flagj   = 1
        flagk   = 0
        var_mtr = jc
      elseif(idir.eq.4) then
        inc     = -param_int(NIJK)
        sens    = -1
        flagi   = 0
        flagj   = 0
        flagk   = 0
        var_mtr = jc
      elseif(idir.eq.5) then
        inc     = param_int(NIJK)*param_int(NIJK+1)
        sens    = 1
        flagi   = 0
        flagj   = 0
        flagk   = 1
        var_mtr = kc
      else
        inc     = -param_int(NIJK)*param_int(NIJK+1)
        sens    = -1
        flagi   = 0
        flagj   = 0
        flagk   = 0
        var_mtr = kc
      endif

      inci =1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      c4 =  5. / 6.
      c5 =  2. / 6.
      c6 = -1. / 6.

      rgp= param_real(CVINF)*(param_real(GAMMA)-1.)  !Cv(gama-1)= R (gas parfait)

      IF (param_int(ITYPZONE).eq.0) THEN

        if(param_int(KFLUDOM).eq.5) then

         do  deep = 0, 1
             !if (deep.eq.1) sens=sens*-1.
         do  f = ideb, ifin

            lc    = facelist( f) + inc*deep
            l1    = lc + inc
            l2    = lc + inc*2
            l0    = lc - inc

            lb    = facelist( f) + inc

            k = (lc-1)/inck -1
            j = (lc-1-(k+1)*inck)/incj -1
            i = (lc-1-(k+1)*inck-(j+1)*incj) -1

            !if (i.lt.1.or.i.gt.param_int(IJKV  )) continue
            !if (j.lt.1.or.j.gt.param_int(IJKV+1)) continue
            !if (k.lt.1.or.k.gt.param_int(IJKV+2)) continue

            lt    = indmtr(i+flagi,j+flagj,k+flagk)
            l     = inddm( i+flagi,j+flagj,k+flagk)

            nm  = l -  inc*sens
            nm2 = l -2*inc*sens
            np  = l +  inc*sens

            tcx  = tijk(lt,1)
            tcy  = tijk(lt,2)
            tcz  = tijk(lt,3)

            surf = sqrt(tcx**2+tcy**2+tcz**2)
            si   = surf

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
            vslp = v1
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm1 = qm
           qp1 = qp
           vslp = v2
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm2 = qm
           qp2 = qp
           vslp = v3
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm3 = qm
           qp3 = qp
           vslp = v4
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm4 = qm
           qp4 = qp
           vslp = v5
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm5 = qm
           qp5 = qp

!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_roe_GD.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_3dfull_i.for"

#include "FastS/Compute/Vit_ent/fludiffer_3dfull_i.for"
          !flu2 = flu2 - tcx*p1p2*0.5
          !flu3 = flu3 - tcy*p1p2*0.5
          !flu4 = flu4 - tcz*p1p2*0.5

       !if(deep.eq.1.and.ndom.eq.0) write(*,*)qp1,qm1,idir
           flu1_old = flu1
           flu2_old = flu2
           flu3_old = flu3
           flu4_old = flu4
           flu5_old = flu5

           drodm(l1,1) =  drodm(l1,1) - flu1*coe(l1,1)*sens
           drodm(l1,2) =  drodm(l1,2) - flu2*coe(l1,1)*sens
           drodm(l1,3) =  drodm(l1,3) - flu3*coe(l1,1)*sens
           drodm(l1,4) =  drodm(l1,4) - flu4*coe(l1,1)*sens
           drodm(l1,5) =  drodm(l1,5) - flu5*coe(l1,1)*sens
           drodm(lc,1) =  drodm(lc,1) + flu1*coe(lc,1)*sens
           drodm(lc,2) =  drodm(lc,2) + flu2*coe(lc,1)*sens
           drodm(lc,3) =  drodm(lc,3) + flu3*coe(lc,1)*sens
           drodm(lc,4) =  drodm(lc,4) + flu4*coe(lc,1)*sens
           drodm(lc,5) =  drodm(lc,5) + flu5*coe(lc,1)*sens

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
            vslp = v1
c#include  "FastS/Compute/Slope/minmod_slope_var.for"
#include  "FastS/Compute/Slope/o1_slope_var.for"
c#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm1 = qm
           qp1 = qp
           vslp = v2
c#include  "FastS/Compute/Slope/minmod_slope_var.for"
#include  "FastS/Compute/Slope/o1_slope_var.for"
c#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm2 = qm
           qp2 = qp
           vslp = v3
c#include  "FastS/Compute/Slope/minmod_slope_var.for"
#include  "FastS/Compute/Slope/o1_slope_var.for"
c#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm3 = qm
           qp3 = qp
           vslp = v4
c#include  "FastS/Compute/Slope/minmod_slope_var.for"
#include  "FastS/Compute/Slope/o1_slope_var.for"
c#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm4 = qm
           qp4 = qp
           vslp = v5
c#include  "FastS/Compute/Slope/minmod_slope_var.for"
#include  "FastS/Compute/Slope/o1_slope_var.for"
c#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm5 = qm
           qp5 = qp

!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_roe_GD.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_3dfull_i.for"

#include "FastS/Compute/Vit_ent/fludiffer_3dfull_i.for"

          if (deep.eq.0) then
            flu1 = 0
            flu2 = tcx*p1p2*0.5
            flu3 = tcy*p1p2*0.5
            flu4 = tcz*p1p2*0.5
            flu5 = 0
          endif

       !if(deep.eq.1) write(*,*)flu1, flu1_old, deep, ndom

          drodm(l1,1) =  drodm(l1,1) + flu1*coe(l1,1)*sens
          drodm(l1,2) =  drodm(l1,2) + flu2*coe(l1,1)*sens
          drodm(l1,3) =  drodm(l1,3) + flu3*coe(l1,1)*sens
          drodm(l1,4) =  drodm(l1,4) + flu4*coe(l1,1)*sens
          drodm(l1,5) =  drodm(l1,5) + flu5*coe(l1,1)*sens
          drodm(lc,1) =  drodm(lc,1) - flu1*coe(lc,1)*sens
          drodm(lc,2) =  drodm(lc,2) - flu2*coe(lc,1)*sens
          drodm(lc,3) =  drodm(lc,3) - flu3*coe(lc,1)*sens
          drodm(lc,4) =  drodm(lc,4) - flu4*coe(lc,1)*sens
          drodm(lc,5) =  drodm(lc,5) - flu5*coe(lc,1)*sens

         enddo
         enddo

       else  !ausm/senseur

         do  f = ideb, ifin

            lc    = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc

            k = (lc-1)/inck -1
            j = (lc-1-(k+1)*inck)/incj -1
            i = (lc-1-(k+1)*inck-(j+1)*incj) -1

            lt    = indmtr(i+flagi,j+flagj,k+flagk)
            l     = inddm( i+flagi,j+flagj,k+flagk)

            nm  = l -  inc*sens
            nm2 = l -2*inc*sens
            np  = l +  inc*sens

            tcx  = tijk(lt,1)
            tcy  = tijk(lt,2)
            tcz  = tijk(lt,3)

            surf = sqrt(tcx**2+tcy**2+tcz**2)
            si   = surf

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
            vslp = v1
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm1 = qm
           qp1 = qp
           vslp = v2
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm2 = qm
           qp2 = qp
           vslp = v3
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm3 = qm
           qp3 = qp
           vslp = v4
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm4 = qm
           qp4 = qp
           vslp = v5
#include  "FastS/Compute/Slope/o3_slope_var.for"
           qm5 = qm
           qp5 = qp

!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_GD.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_3dfull_i.for"

         ! modification de vitesse normale par ajout
         ! de stabilisation de type Rhie-Chow
          c   = rgp*gam*rop(l +v5)  !c^2
          son =sqrt(qn1*qn1 / c)
          tam =c3*son+surf
          tam1=max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
          u   =0.25*(qn1+qn2)-tam1*(p2-p1)
          tdu = max(abs(u),c1*surf)

          !Calcul du flux total
          p1p2= (p1+p2)*0.5*0.

#include "FastS/Compute/Vit_ent/fluvector_3dfull_i.for"

           drodm(l1,1) =  drodm(l1,1) - flu1*coe(l1,1)*sens
           drodm(l1,2) =  drodm(l1,2) - flu2*coe(l1,1)*sens
           drodm(l1,3) =  drodm(l1,3) - flu3*coe(l1,1)*sens
           drodm(l1,4) =  drodm(l1,4) - flu4*coe(l1,1)*sens
           drodm(l1,5) =  drodm(l1,5) - flu5*coe(l1,1)*sens
         enddo
       endif

      ELSEIF  (param_int(ITYPZONE).eq.2) THEN !3Dcart

         do  f = ideb, ifin

            l     = facelist( f)
            l1    = facelist( f) + inc
            l2    = facelist( f) + inc*2
            l0    = facelist( f) - inc
            lt    = 1
            
            tcx  = tijk(lt,ic)*ci_mtr
            tcy  = tijk(lt,jc)*cj_mtr
            tcz  = tijk(lt,kc)*ck_mtr

            !u = rop(l , 2)*tcx + rop(l , 3)*tcy + rop(l , 4)*tcz
            !u1= rop(l1, 2)*tcx + rop(l1, 3)*tcy + rop(l1, 4)*tcz
            !u2= rop(l2, 2)*tcx + rop(l2, 3)*tcy + rop(l2, 4)*tcz
            !u0= rop(l0, 2)*tcx + rop(l0, 3)*tcy + rop(l0, 4)*tcz

            !flux(1) = flux(1) +  u*sens*rop(l, 1)
            !flux(2) = flux(2) + u1*sens*rop(l1,1)
            !flux(3) = flux(3) + u2*sens*rop(l2,1)
            !flux(7) = flux(7) + u0*sens*rop(l0,1)
            !!!flux(1) = flux(1) +  u*sens
            !!!flux(2) = flux(2) + u1*sens
            !!!flux(3) = flux(3) + u2*sens
            !!!flux(7) = flux(7) + u0*sens
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

            !u = c5*rop(l, 2) + c4*rop(l1,2) + c6*rop(l2,2)
            !v = c5*rop(l, 3) + c4*rop(l1,3) + c6*rop(l2,3)

            k = l/inck
            j = (l-k*inck)/incj
            i = l-k*inck - j*incj
            qn = (u*tcx+v*tcy)

            !flux(1)=flux(1)+(rop(l ,2)*tcx+rop(l ,3)*tcy)*sens*rop(l,1)
            !flux(2)=flux(2)+(rop(l1,2)*tcx+rop(l1,3)*tcy)*sens*rop(l1,1)
            !flux(3)=flux(3)+(rop(l2,2)*tcx+rop(l2,3)*tcy)*sens*rop(l2,1)
            !flux(7)=flux(7)+(rop(l0,2)*tcx+rop(l0,3)*tcy)*sens*rop(l0,1)
            !flux(4) = flux(4) +  tcx + tcy

            !flux(5) = flux(5) +  (tcx + tcy + tcz)*sens

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
c            flux(6) = flux(6) +  tj(lt,2)*rop(l,1)*rop(l,5)*rgp0.000512834373688 -0.000105837563190
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
