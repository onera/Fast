c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 39 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bflwallmodel(ndom, idir, iptparam, nitcfg, neq_mtr,
     &                   param_int, param_real, incijk, ind_loop,
     &                   rop, drodm, tijk, ventijk, xmut, sampling)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  bflco0k.sf/pechier1/e1/i1
c
c_DC  DATE_C : Octobre 1990 -- GUILLEN / DORMIEUX
c
c     HISTORIQUE
c_H    Iteration     1 ___ 1997-06-05 15:34:15  pechier
c_H    Issu de bflco0k.sf/lepape1/e1/i4
c
c_U   USER : GUILLEN
c
c     ACT
c_A    Calcul de flux frontieres aux interfaces solides du domaine,
c_A    la pression est obtenue par relations de compatibilite.
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    processeur domaine
c
c     INP
c_I    ndom   : numero du domaine
c_I    idir   : direction du traitement
c_I    ndimdx : dimension 3D
c_I    tijk   : tableaux des normales
c_I    flu    : tableau des flux explicites
c_I    q      : variables d etats aux interfaces
c
c     OUT
c
c     I/O
c_/    flu    : increment des variables conservatives
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,neq_mtr,incijk, idir, ind_loop(6), param_int(0:*)

      REAL_E   rop( param_int(NDIMDX)     * param_int(NEQ)     )
      REAL_E drodm( param_int(NDIMDX)     * param_int(NEQ)     )
      REAL_E ventijk( param_int(NDIMDX_VENT)* param_int(NEQ_VENT))

      REAL_E  tijk( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) )

      REAL_E xmut(*)

      REAL_E param_real(0:*)

C var loc
      INTEGER_E im,jm,km,ijkm,l,l0,iadrf,m,i,j,k,lj,ic,jc,kc,
     & kc_vent,v1,v2,v3,v4,v5,v6,vmtr,vven,shift,iptparam,
     & lt,lven,lij,ltij,lvij,iter,sampling,input_value,flag,
     & exchange, nitcfg, ilaw

      REAL_E p,r,u,v,w,qen,ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,
     & ck_mtr_vent, u_int,sens,rgp,flu1,flu2,flu3,flu4,flu5,flu6,
     & tcx,tcy,tcz,dist,utau,unorm,surf,aa,yplus,l1,l2,l3,f,tp,fp,tauw,
     & ut,vt,wt,Twall,uplus,pr,prtur,ratio,dt,pf,Tplus,qwall,cp,cond
 
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

      rgp    = param_real( CVINF )*(param_real( GAMMA )-1.)  !Cv(gama-1)= R (gas parfait)

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v6 = 5*param_int(NDIMDX)

      vmtr  =   param_int(NDIMDX_MTR)
      vven  =   param_int(NDIMDX_VENT)

!     sampling = iptparam
!     sampling = 0
!     if (flag.eq.1) sampling = -input_value !add Kawai and Larsson values
      exchange = 1 - sampling  ! index correction for identifying the correct cell

      sens  = 1.
      shift = 0
      if(mod(idir,2).eq.0) then
        shift = 1 
        exchange = 1 - exchange
      endif

      shift =incijk*shift
      exchange =incijk*exchange

      !!  a mettre dans Fast
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      !correction monoindice
      ic      = ic -1
      jc      = jc -1
      kc      = kc -1
      kc_vent = kc_vent -1

      IF(param_int(NEQ).eq.5) THEN

         DO k = ind_loop(5), ind_loop(6)
         DO j = ind_loop(3), ind_loop(4)
         DO i = ind_loop(1), ind_loop(2)

          l     = inddm(  i, j, k)
          lt    = indmtr( i, j, k)
          lven  = indven( i, j, k)

          iadrf = l  - incijk
          l0    = l  - shift
          m     = l  - exchange

          tcx = tijk(lt +vmtr*ic)*ci_mtr
          tcy = tijk(lt +vmtr*jc)*cj_mtr
          tcz = tijk(lt +vmtr*kc)*ck_mtr

          surf = sqrt(tcx*tcx+tcy*tcy+tcz*tcz)

          ! Calculate tangential velocity and its norm
          u = rop(m+v2)
          v = rop(m+v3)
          w = rop(m+v4)

          u_int = tcx*u +tcy*v +tcz*w ! normal velocity

          ut = u-u_int*tcx
          vt = v-u_int*tcy
          wt = w-u_int*tcz

          unorm = sqrt(ut*ut+vt*vt+wt*wt)

          r     = 0.5*(rop(l+v1)+rop(iadrf+v1))

          p = 0.5*(rop(l+v5)*rop(l+v1)+rop(iadrf+v5)*rop(iadrf+v1))*rgp

#include  "FastS/BC/law_velocity.for"
 
!           if (prtur.gt.0.1) then
!              Twall = 0.5*(rop(l+v5)+rop(iadrf+v5))
!#include  "FastS/BC/law_temperature.for"
!           end if

          flu1 = 0.
          flu2 = tcx*p - (tauw*rop(m+v2)/unorm)*surf
          flu3 = tcy*p - (tauw*rop(m+v3)/unorm)*surf
          flu4 = tcz*p - (tauw*rop(m+v4)/unorm)*surf
          flu5 = qwall*surf

#include  "FastS/Compute/assemble_drodm_corr.for"
       enddo
       enddo
       enddo

      ENDIF

      end
