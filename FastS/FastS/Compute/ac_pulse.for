c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine ac_pulse(ndom, param_int, ind_loop,
     &                    temps, staref, 
     &                    drodm,rop,x,y,z,vol,
     &                    x0,y0,z0,coefa,coefb,amp,per,phi)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  ac_pulse.sf
c
c_DC  DATE_C : Janvier 2002 - M. Terracol 
c
c     HISTORIQUE
c
c     ACT
c_A    Generation d''un pulse acoustique 
c
c     VAL
c


c     INP
c_I    ndom   : numero du domaine
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, ind_loop(6),param_int(0:*)


      REAL_E   rop( param_int(NDIMDX)   , param_int(NEQ) )
      REAL_E drodm( param_int(NDIMDX)   , param_int(NEQ) )
      REAL_E   vol( param_int(NDIMDX_MTR) )
      REAL_E     x( param_int(NDIMDX_XYZ) )
      REAL_E     y( param_int(NDIMDX_XYZ) )
      REAL_E     z( param_int(NDIMDX_XYZ) )

      REAL_E temps,staref(*), x0,y0,z0,coefa,coefb,amp,per,phi
      
C
C Var loc
      INTEGER_E i, j, k,l,lij,inci,incj,inck,ios,lt,lx,ltij,lxij, lvo

      REAL_E pi,gam,alpha,omega,omt,eps,xx,yy,zz,r2,r,u,v,w,q,p,sro,
     & sroe,roe,seuil,gamma,rgp,gam2,gam3,diff

#include "FastS/formule_param.h"
#include "FastS/formule_xyz_param.h"
#include "FastS/formule_mtr_param.h"

      inci = 1
      incj = param_int(NIJK_XYZ)
      inck = param_int(NIJK_XYZ)*param_int(NIJK_XYZ+1)


      pi     = 4.*atan(1.)
      gamma  = staref(1)
      rgp    = staref(2)*(gamma-1.)  !Cv(gama-1)= R (gas parfait)
      gam2   =  gamma*rgp
      gam3   =  gam2/(gamma-1.)

      alpha  = log(coefa)/(coefb*coefb)
      omega  = 2.*pi/per
      omt    = omega*temps+phi
      eps    = amp*sin(omt)

      seuil = exp(-alpha*30.)
 
      !on split pour vectoriser en compilo v11
      IF(param_int(ITYPZONE).eq.2) THEN  !3D cart


        do 10 k = ind_loop(5), ind_loop(6)
        do 10 j = ind_loop(3), ind_loop(4)

          lij  = inddm(  ind_loop(1) , j, k)
          lxij = lij -  indcg( ind_loop(1) , j, k)
          lt   = indmtr( 1 , 1, 1)
          lvo  = lt

!DEC$ IVDEP
          do 10 l = lij,  lij + ind_loop(2) - ind_loop(1)

            lx = l  - lxij


             xx = 0.125*( x(lx)      +x(lx+inci)     +x(lx+inci+incj)  
     &                + x(lx+incj) 
     &                + x(lx+inck) +x(lx+inci+inck)+x(lx+inci+incj+inck)
     &                + x(lx+incj+inck) )

             yy =0.125*(y(lx)+y(lx+inci)+y(lx+inci+incj)+y(lx+incj)
     &               +y(lx+inck)  +y(lx+inci+inck)+y(lx+inci+incj+inck)
     &               +y(lx+incj+inck))
             zz = z0

            r2 = (xx-x0)*(xx-x0) + (yy-y0)*(yy-y0) + (zz-z0)*(zz-z0)

c            r2=r2*r2
c            sro = eps*vol(lvo)/(1.+r2)
            if(r2.ge.30) then
              sro = eps*seuil*vol(lvo)
            else
             sro = eps*exp(-alpha*r2)*vol(lvo)
            endif
        
            !c0  = gam2*rop(l,5)
            !sp  = c0*sro
            sroe=  rop(l,5)*sro*gam3
  
            drodm(l,1) = drodm(l,1) + sro
            drodm(l,5) = drodm(l,5) + sroe

   10 continue


      ELSE


      do 50 k = ind_loop(5), ind_loop(6)
      do 50 j = ind_loop(3), ind_loop(4)

        lij  = inddm(  ind_loop(1) , j, k)
        lxij = lij - indcg( ind_loop(1) , j, k)
        ltij = lij - indmtr(ind_loop(1) , j, k)

!DEC$ IVDEP
        do 50 l = lij,  lij + ind_loop(2) - ind_loop(1)

          lx = l  - lxij
          lt = l  - ltij
          lvo= lt

        xx = 0.125*( x(lx)      + x(lx+inci)     + x(lx+inci+incj)     
     &             + x(lx+incj) 
     &             + x(lx+inck) + x(lx+inci+inck)+ x(lx+inci+incj+inck)
     &             + x(lx+incj+inck) )

        yy =0.125*(y(lx)+y(lx+inci)+y(lx+inci+incj)+y(lx+incj)
     &            +y(lx+inck)      +y(lx+inci+inck)+y(lx+inci+incj+inck)
     &            +y(lx+incj+inck))

        zz = z0


c        xx = -49.5 + i
c        yy = -49.5 + j
c        zz = -49.5 + k
c
        r2 = (xx-x0)*(xx-x0) + (yy-y0)*(yy-y0) + (zz-z0)*(zz-z0)

c        r2=r2*r2
c        sro = eps*vol(lvo)/(1.+r2)
        if(r2.ge.30) then
        sro = eps*seuil*vol(lvo)
        else
        sro = eps*exp(-alpha*r2)*vol(lvo)
        endif
        
        !c0  = gam2*rop(l,5)
        !sp  = c0*sro
        sroe=  rop(l,5)*sro*gam3

        !write(*,'(a,f18.12,3i7)')'puls',vol(lvo),lt,j,k

        drodm(l,1) = drodm(l,1) + sro
        drodm(l,5) = drodm(l,5) + sroe

   50 continue

      ENDIF

      end
