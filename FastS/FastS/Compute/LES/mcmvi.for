c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine mcmvi(ndom, param_int, param_real, neq_rot,
     &                 ind_loop,
     &                 xmut, rop, rot)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  mcmvi.sf
c
c_DC  DATE_C : Fevrier 1999 -- DAUX
c     HISTORIQUE
c     ACT
c_A    Calcul de la viscosite de sous-maille, modele MCM
c
c     VAL
c
c     INP
c_I    ndom   : numero du domaine
c_I    ro     : variables conservatives 
c
c     OUT
c_O    xmules : viscosite turbulente de sous-maille
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"
 
      INTEGER_E ndom, neq_rot, ind_loop(6), param_int(0:*)

      REAL_E xmut( param_int(NDIMDX) )
      REAL_E rot ( param_int(NDIMDX) * neq_rot)
      REAL_E rop( param_int(NDIMDX)  * param_int(NEQ) )

      REAL_E param_real(0:*)

C var loc
      INTEGER_E i,j,k,l,lt, lvo, ltij,lij, inci,incj,inck,
     & l0,l1,l2,l3, l4,l5,l6,l7,l8,var,v1,v2,v3,v4,v5

      REAL_E pi, angle,cste,vref,vort,vortmoy,vortflu,r1,r2,fs,xkc_l,fv,
     & c111,c110,c100,c000,c1,c0,temp01,cmus1,coesut,t1

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

      v1 =   0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      cmus1  =    param_real( CS )
      temp01 = 1./param_real( TEMP0 )
      coesut =    param_real( XMUL0 ) * (1.+cmus1*temp01)*sqrt(temp01)

      pi    = 4. * atan(1.)
      angle = 45 * pi / 360.
      vref  = 1. / max(tan(angle)**2, 1.e-20)
      cste  = 0.08

c     c1 = 0.25
c     c0 = 0.5
c     c1 = 1./6.
c     c0 = 4./6.
c     c1 = 1./3.
c     c0 = 1./3.

      !!!Funk Terracol value
c      rap= sqrt(6.)
c      c1 = rap**2/24.
c      c0 = (12.-rap**2)/12.
      c1 = 0.25
      c0 = 0.5

      c111 = c1*c1*c1
      c110 = c1*c1*c0
      c100 = c1*c0*c0
      c000 = c0*c0*c0

      !angle = mcmangle(ndom) * pi / 180.
      !angle = mcmangle(ndom) * pi / 360.

      DO k = ind_loop(5), ind_loop(6)
       DO j = ind_loop(3), ind_loop(4)

#include "FastS/Compute/loopI_begin.for"


             !Srate du funk calculer dans cp_rot_sij.for
             xmut(l) =   xmut(l)**0.25

             !fitrage VelocityX
             var = v2
#include "FastS/Compute/LES/filtre_vitesse.for"
             xkc_l = (fv - rop(l +var) )*(fv - rop(l +var) )
c             if(k.eq.20.and.j.le.2.and.l-lij.le.3) then
c        write(*,'(a,2f19.15,3i4)')'vi',
c     &           fv , rop(l +var)  ,l-lij,j,k
c             endif

             !xmut(l)=fv

             !fitrage VelocityY
             var = v3
#include "FastS/Compute/LES/filtre_vitesse.for"
             xkc_l = xkc_l + (fv - rop(l +var) )*(fv - rop(l +var) )

c             if(k.eq.20.and.j.le.2.and.l-lij.le.3) then
c        write(*,'(a,2f19.15,3i4)')'vj',
c     &           fv , rop(l +var)  ,l-lij,j,k
c             endif
c             rop(l)=fv

             !fitrage VelocityZ
             var = v4
#include "FastS/Compute/LES/filtre_vitesse.for"
             xkc_l = xkc_l + (fv - rop(l +var) )*(fv - rop(l +var) )

c             if(k.eq.20.and.j.le.2.and.l-lij.le.3) then
c        write(*,'(a,2f19.15,3i4)')'vk',
c     &           fv , rop(l +var), l-lij,j,k
c        write(*,'(a,3f19.15)')'mu',          xmut(l)
c             endif
c             rop(l+v5)=fv


             xkc_l     = 1./sqrt(2.) * sqrt ( xkc_l)

             xmut(l) = cste * rop(l +v1) * sqrt(xkc_l) * xmut(l)


             !filtrage rotX
             var = v1
#include "FastS/Compute/LES/filtre_rot.for"
             vort    = rot(l +var)*rot(l +var)
             vortmoy = fv*fv
             vortflu = (fv - rot(l +var) )*(fv - rot(l +var) )

             !filtrage rotY
             var = v2
#include "FastS/Compute/LES/filtre_rot.for"
             vort    = vort    + rot(l +var)*rot(l +var)
             vortmoy = vortmoy + fv*fv
             vortflu = vortflu + (fv - rot(l +var))*(fv - rot(l +var))

             !filtrage rotZ
             var = v3
#include "FastS/Compute/LES/filtre_rot.for"
             vort    = vort    + rot(l +var)*rot(l +var)
             vortmoy = vortmoy + fv*fv
             vortflu = vortflu + (fv - rot(l +var))*(fv - rot(l +var))


             !selection function
             vort    = sqrt ( vort)
             vortmoy = sqrt ( vortmoy)
             vortflu = sqrt ( vortflu) 

!     tan^2(teta/2)=r1/r2
             r1        =     -(vort - vortmoy)**2 + vortflu**2
             r2        = max( (vort + vortmoy)**2 - vortflu**2, 1.e-20 )
             fs        = 1.65*min( 1., (vref*r1/r2)**2 )
 


             ! mut = mulam + mu_sgs
             t1     = rop(l + v5)
             xmut(l) = fs * xmut(l) + coesut * sqrt(t1)*t1/(t1+cmus1)
             !xmut(l) =  xmut(l) + coesut * sqrt(t1)*t1/(t1+cmus1)

c             if(k.le.3.and.k.ge.1.and.j.eq.81.and.l-lij.eq.10) then
c        write(*,'(a,5f19.12,3i4)')'vi',
c     &             xmut(l),fs,vort,vortmoy,vortflu,l-lij,j,k
c     &             xmut(l),fs,r1,r2,vortflu,l-lij,j,k
c             endif

          enddo
       ENDDO !do j
      ENDDO !do k

      end
