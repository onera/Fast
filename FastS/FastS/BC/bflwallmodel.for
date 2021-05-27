c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 39 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bflwallmodel(ndom, ithread,idir,nbdata,nitcfg, neq_mtr,
     &                   data, 
     &                   param_int, param_real, incijk, ind_loop,
     &                   rop, drodm, tijk, ventijk, xmut,x,y,z)
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

      INTEGER_E ndom,neq_mtr,incijk, idir, nbdata,ithread,
     & ind_loop(6), param_int(0:*)

      REAL_E   rop( param_int(NDIMDX)     * param_int(NEQ)     )
      REAL_E drodm( param_int(NDIMDX)     * param_int(NEQ)     )
      REAL_E ventijk( param_int(NDIMDX_VENT)* param_int(NEQ_VENT))

      REAL_E  tijk( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) )

      REAL_E xmut(*),x(*),y(*),z(*)

      REAL_E param_real(0:*)
      REAL_E data(0:*)

C var loc
      INTEGER_E im,jm,km,ijkm,l,l0,iadrf,m,i,j,k,lj,ic,jc,kc,
     & kc_vent,v1,v2,v3,v4,v5,v6,vmtr,vven,shift,shiftvent,
     & lt,lven,lij,ltij,lvij,iter,sampling,
     & exchange, nitcfg, ilaw,mvent,lx,lr

      REAL_E p,r,u,v,w,qen,ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,
     & ck_mtr_vent, u_int,sens,rgp,flu1,flu2,flu3,flu4,flu5,flu6,
     & nx, ny,nz,utau2,seuil,mobile_coef,ventx,venty,ventz,unorm1,
     & tcx,tcy,tcz,dist,utau,unorm,surf,aa,yplus,l1,l2,l3,f,tp,fp,tauw,
     & ut,vt,wt,Twall,uplus,pr,prtur,ratio,dt,pf,Tplus,qwall,cp,cond,
     & wall_ut,wall_vt,wall_wt,wall_int,co,si,rot1,rot2,ra,un,vn,wn,umod
 
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_xyz_param.h"
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

      sampling = param_int(WM_SAMPLING)

      exchange = 1 - sampling  ! index correction for identifying the correct cell

      sens  = 1.
      shift = 0
      if(mod(idir,2).eq.0) then
        sens  = -1.
        shift = 1 
        exchange = 1 - exchange
      endif

      !write(*,*)'shiftvent', idir/2,NIJK_VENT

      shift =incijk*shift
      shiftvent = param_int(NIJK_VENT + idir/2)* exchange


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

      mobile_coef = 1.
      if(nbdata.ne.0) mobile_coef = data(0)

      !write(*,*)'loop', ind_loop, idir

      IF(param_int(NEQ).eq.5) THEN

        ilaw = param_int(WM_FUNCTION)

        if (ilaw.eq.0) then  !musker

           DO k = ind_loop(5), ind_loop(6)
           DO j = ind_loop(3), ind_loop(4)
           DO i = ind_loop(1), ind_loop(2)

            l     = inddm(  i, j, k)
            lt    = indmtr( i, j, k)
            lven  = indven( i, j, k)

            iadrf = l  - incijk
            l0    = l  - shift
            m     = l  - exchange
            mvent = lven  - shiftvent

            tcx = tijk(lt +vmtr*ic)*ci_mtr
            tcy = tijk(lt +vmtr*jc)*cj_mtr
            tcz = tijk(lt +vmtr*kc)*ck_mtr

            surf = sqrt(tcx*tcx+tcy*tcy+tcz*tcz)
            surf = max(surf,1e-30)

            nx = tcx/surf
            ny = tcy/surf
            nz = tcz/surf

            ventx      = ventijk(lven              )*c_ale
            venty      = ventijk(lven +vven        )*c_ale
            ventz      = ventijk(lven +vven*kc_vent)*ck_vent*c_ale

          qen =(  ventijk(lven              )*tcx
     &           +ventijk(lven +vven        )*tcy
     &           +ventijk(lven +vven*kc_vent)*tcz*ck_vent
     &         )*c_ale

            ! Calculate tangential velocity and its norm
            u = rop(m+v2) 
            v = rop(m+v3) 
            w = rop(m+v4) 

            u_int = nx*u + ny*v + nz*w ! normal velocity

            un = u_int*nx
            vn = u_int*ny
            wn = u_int*nz

            wall_int = nx*ventx + ny*venty + nz*ventz

            wall_ut = ventx-wall_int*nx
            wall_vt = venty-wall_int*ny
            wall_wt = ventz-wall_int*nz

            !Modif Ivan
            ut = u-u_int*nx - wall_ut
            vt = v-u_int*ny - wall_vt
            wt = w-u_int*nz - wall_wt

            unorm = sqrt(ut*ut+vt*vt+wt*wt)
c      if(ndom.eq.3.and.j.eq.1.and.i.eq.150.and.k.eq.30)
c     & write(*,'(a,4f20.15)')'U',ut,u,u_int, wall_ut

            r     = 0.5*(rop(l+v1)+rop(iadrf+v1))

            p =0.5*(rop(l+v5)*rop(l+v1)+rop(iadrf+v5)*rop(iadrf+v1))*rgp

            dist = xmut(m+param_int(NDIMDX)) !wall distance

            ! Wall function for the velocity proposed in
            ! "Explicit expression for the smooth wall velocity
            ! distribution in a turbulent boundary layer"
            ! AIAA Journal. Musker (1979)

            ! Initial prediction based on linear profile (u+=y+)
            if(sampling.eq.1) then
              utau = sqrt( (xmut(l)*unorm)/(rop(m+v1)*dist) ) 
            else
              utau=((unorm**7)*xmut(l)/(rop(m+v1)*dist)/(8.3**7))**0.125
            endif

c      if(ndom.eq.3.and.j.eq.1.and.i.eq.150.and.k.eq.30)
c     & write(*,'(a,6f20.15,i4)')'utau_init',r*utau**2,unorm,
c     & rop(m+v1),rop(m+v2),dist,iter

            ! Newton method to determine utau using Musker's wall law
            iter = 0 
            !aa    = rop(m+v1)*dist/xmut(m)
            aa    = rop(m+v1)*dist/xmut(l)
            yplus = aa*utau
            l1    = (yplus+10.6)**9.6
            l2    = yplus*(yplus-8.15) + 86
            l3    = (2*yplus-8.15)/16.7
            f = 5.424*atan(l3) + log10(l1/(l2*l2)) - unorm/utau - 3.52
c      if(ndom.eq.3.and.j.eq.1.and.i.eq.150.and.k.eq.30)
c     & write(*,'(a,4f20.15)')'f init',f,utau,yplus,- unorm/utau - 3.52

            do while (abs(f).gt.1e-6.and.iter.le.50)
 
             iter = iter + 1
 
             !aa    = rop(m+v1)*dist/xmut(m)
             aa    = rop(m+v1)*dist/xmut(l)
             yplus = aa*utau
             l1    = (yplus+10.6)**9.6
             l2    = yplus*(yplus-8.15) + 86
             l3    = (2*yplus-8.15)/16.7
 
             tp =aa*(9.6*((yplus + 10.6)**8.6)- l1/l2*( 4*yplus-16.30) )
 
             fp = 0.649580838323*aa/(1 + l3*l3 )  +  tp/(l1*2.302585093)
     &        + unorm/(utau*utau)
 
             !blinder fp
             utau = abs(utau - f/fp)

             !aa    = rop(m+v1)*dist/xmut(m)
             aa    = rop(m+v1)*dist/xmut(l)
             yplus = aa*utau
             l1    = (yplus+10.6)**9.6
             l2    = yplus*(yplus-8.15) + 86
             l3    = (2*yplus-8.15)/16.7

             f = 5.424*atan(l3) + log10(l1/(l2*l2)) - unorm/utau - 3.52

            end do

            !if(iter.eq.51)write(*,*)'Newton musker', iter,ndom

             !if(ndom.eq.3.and.j.eq.1.and.i.eq.150.and.k.eq.30)
c             if(ndom.eq.0.and.j.eq.1.and.i.eq.10.and.k.eq.30)
c     &       write(*,'(a,5f20.15,2i4)')'utau_finale',r*utau**2,unorm,
c     &       f,rop(m+v2),dist,iter,k


            utau2 = utau**2
            tauw = r*utau2

c      if(k.eq.10.and.j.le.15.and.idir.eq.1)
c     & write(*,'(a,f20.10,i5)')'verI', sqrt(utau2),j
c      if(k.eq.10.and.i.le.5.and.idir.eq.3)
c     & write(*,'(a,f20.10,i5)')'verJ', sqrt(utau2),i

            !unorm = sqrt(rop(l0+v2)**2+rop(l0+v3)**2+rop(l0+v4)**2)
            !unorm = max(1.,unorm)

            u     = 0.5*(rop(l+v2)+rop(iadrf+v2))
            v     = 0.5*(rop(l+v3)+rop(iadrf+v3))
            w     = 0.5*(rop(l+v4)+rop(iadrf+v4))

c          lx =indcg(i, j, k)
c
c          ra  = sqrt( y(lx)*y(lx) + z(lx)*z(lx))
c          ra  = max(ra, 0.00000001)
c          co = z(lx)/ra
c          si = y(lx)/ra
c          rot1 = (-venty*co+ventz*si)
c          rot2 = (-v*co+w*si)
c          if(idir.eq.3) write(*,'(4f20.12,2i4)')
c     &   rot1, rot2, ra*3.1415*11543./30.,ra,j,k
c          if(idir.eq.3) write(*,'(a,3f20.12,2i4)')'JJ',
c     &   u_int/surf, (tcx*u +tcy*v +tcz*w)/surf, qen/surf,i,k
c          if(idir.eq.1) write(*,'(a,3f20.12,2i4)')'II',
c     &   u_int/surf, (tcx*u +tcy*v +tcz*w)/surf, qen/surf,j,k


            !determination vitesse normale interface
            u_int= tcx*u +tcy*v +tcz*w -qen 


c     & - (tauw*ut/unorm)*surf*sens,tauw,ut/unorm,surf,dist

            flu1= 0.
            flu2= tcx*p   - (tauw*ut/unorm)*surf*sens
            flu3= tcy*p   - (tauw*vt/unorm)*surf*sens
            flu4= tcz*p   - (tauw*wt/unorm)*surf*sens
            flu5= 0.
            !flu1= 0.
            !flu2= tcx*p  + u_int*u*r - (tauw*ut/unorm)*surf*sens
            !flu3= tcy*p  + u_int*v*r - (tauw*vt/unorm)*surf*sens
            !flu4= tcz*p  + u_int*w*r - (tauw*wt/unorm)*surf*sens
            !flu5= p*qen              - qwall*surf*sens

            yplus = rop(l+v1)*utau*xmut(l+v2)/xmut(l)
            umod = utau*(5.424*atan((2.*yplus - 8.15)/16.7)
     &            + log10((yplus + 10.6)**9.6
     &            /(yplus**2 - 8.15*yplus + 86.)**2) - 3.52)

            umod = abs(umod)
            lr   = l + incijk

cc            if(nitcfg.eq.1) then
c             rop(l+v1) = rop(lr+v1)
c             rop(l+v5) = rop(lr+v5)
cc            endif
c            rop(l+v2) = umod*ut/unorm + xmut(l+v2)/xmut(m+v2)*un
c            rop(l+v3) = umod*vt/unorm + xmut(l+v2)/xmut(m+v2)*vn
c            rop(l+v4) = umod*wt/unorm + xmut(l+v2)/xmut(m+v2)*wn


c           rop(l+v2)     = tauw*ut/unorm*xmut(l+v2)/xmut(l)
c           rop(iadrf+v2) = -rop(l+v2)

c      if(ndom.eq.1.and.j.eq.1.and.i.eq.130.and.k.eq.30)
c     & write(*,'(a,6f20.15,i6)')'wall',
c     & flu2,rop(l+v2),-tauw*ut/unorm,rop(m+v2),
c     & rop(l+v2),tauw,xmut(l),xmut(l+v2) ,rop(m+v2), 
c     & rop(l+v4)

#include  "FastS/Compute/assemble_drodm_corr.for"
           enddo
           enddo
           enddo

        elseif(ilaw.eq.1) then !powerlaw
 
           DO k = ind_loop(5), ind_loop(6)
           DO j = ind_loop(3), ind_loop(4)
           DO i = ind_loop(1), ind_loop(2)

            l     = inddm(  i, j, k)
            lt    = indmtr( i, j, k)
            lven  = indven( i, j, k)

            iadrf = l  - incijk
            l0    = l  - shift   !1er point au dessus paroi
            m     = l  - exchange
            !mvent = lven  - shiftvent

            tcx = tijk(lt +vmtr*ic)*ci_mtr
            tcy = tijk(lt +vmtr*jc)*cj_mtr
            tcz = tijk(lt +vmtr*kc)*ck_mtr

            surf = sqrt(tcx*tcx+tcy*tcy+tcz*tcz)
            surf = max(surf,1e-30)

            nx = tcx/surf
            ny = tcy/surf
            nz = tcz/surf

            ventx      = ventijk(lven              )*c_ale
            venty      = ventijk(lven +vven        )*c_ale
            ventz      = ventijk(lven +vven*kc_vent)*ck_vent*c_ale


          qen =(  ventijk(lven              )*tcx
     &           +ventijk(lven +vven        )*tcy
     &           +ventijk(lven +vven*kc_vent)*tcz*ck_vent
     &         )*c_ale

            ! Calculate tangential velocity and its norm
            u = rop(m+v2) 
            v = rop(m+v3) 
            w = rop(m+v4) 

            u_int = nx*u + ny*v + nz*w ! normal velocity

            wall_int = nx*ventx + ny*venty + nz*ventz
            wall_ut = ventx-wall_int*nx
            wall_vt = venty-wall_int*ny
            wall_wt = ventz-wall_int*nz

            !Modif Ivan
            ut = u-u_int*nx - wall_ut
            vt = v-u_int*ny - wall_vt
            wt = w-u_int*nz - wall_wt

            unorm = sqrt(ut*ut+vt*vt+wt*wt)

            r     = 0.5*(rop(l+v1)+rop(iadrf+v1))

            p =0.5*(rop(l+v5)*rop(l+v1)+rop(iadrf+v5)*rop(iadrf+v1))*rgp

            dist = xmut(m+param_int(NDIMDX)) !wall distance

            ! Power law proposed in
            ! "Large-eddy simulation of turbulent flow over and around
            ! a cube in a plate channel"
            ! Werner & Wengle (1991)
      
            !u+=y+  si y+ < 11.8    utau = sqrt (Mu * Unorm/ density/ dist)
            !u+=8.3 (y+)**1/7  sinon  utau =( Unorm**7/8.3**7 * mu/density/dist)**1/8
    
            seuil = 8.3**(7./3.)*xmut(m)/(rop(m+v1)*dist)

c             if(ndom.eq.4.and.idir.eq.1)
c     &   write(*,'(a,3f15.10,6i7)')'seuil',unorm,xmut(m),
c     &     rop(m+v1), dist,i,j,k,lven,vven
c     &   write(*,'(a,5f15.10,4i4)')'seuil',seuil,unorm,xmut(m),
c     &     rop(m+v1), dist,i,j,k,ndom

            if(unorm.le.seuil) then
              utau2 = xmut(m)/rop(m+v1)*unorm/dist
            else
              utau2=( (unorm**7)*xmut(m)/(rop(m+v1)*dist)/(8.3**7)
     &              )**0.25
            endif

            tauw = r*utau2

            qwall = 0 ! adiabatic wall
 
            u     = 0.5*(rop(l+v2)+rop(iadrf+v2))
            v     = 0.5*(rop(l+v3)+rop(iadrf+v3))
            w     = 0.5*(rop(l+v4)+rop(iadrf+v4))


            !determination vitesse normale interface
            u_int= tcx*u +tcy*v +tcz*w -qen

            unorm = sqrt(rop(l0+v2)**2+rop(l0+v3)**2+rop(l0+v4)**2)
            unorm = max(1.,unorm)

c      if(k.eq.10.and.j.le.15.and.idir.eq.1)
c     & write(*,'(a,f20.10,i5)')'verI', sqrt(utau2),j
c      if(k.eq.10.and.i.le.5.and.idir.eq.3)
c     & write(*,'(a,f20.10,i5)')'verJ', sqrt(utau2),i
c      if(k.eq.10.and.j.eq.10)write(*,'(a,9f18.14)')'verI',rop(m+v2),
c     & rop(l+v2),u_int,qen/surf,venty,ventz,nx,ny,nz,idir
c      if(k.eq.10.and.j.eq.10)write(*,'(a,7f19.14)')'verI',
c     & u_int/surf,qen/surf,venty,ventz,nx,ny,nz,idir

c      if(k.eq.10.and.i.eq.10)write(*,'(a,7f19.14)')'verJ',
c     & u_int/surf,qen/surf,venty,ventz,nx,ny,nz

c      if(k.eq.10.and.j.eq.10.and.idir=1)write(*,*)'verifI',rop(m+v2),
c     & rop(l+v2),m,l,idir
c      if(idir.eq.4)write(*,*)'verifJ',rop(m+v2),

            flu1= 0.
            !flu2= tcx*p  + u_int*u*r - (tauw*rop(l0+v2)/unorm)*surf*sens
            !flu3= tcy*p  + u_int*v*r - (tauw*rop(l0+v3)/unorm)*surf*sens
            !flu4= tcz*p  + u_int*w*r - (tauw*rop(l0+v4)/unorm)*surf*sens
            flu2= tcx*p  + u_int*u*r - (tauw*ut/unorm)*surf*sens
            flu3= tcy*p  + u_int*v*r - (tauw*vt/unorm)*surf*sens
            flu4= tcz*p  + u_int*w*r - (tauw*wt/unorm)*surf*sens
            flu5= p*qen              - qwall*surf*sens
            !flu5= 0.


#include  "FastS/Compute/assemble_drodm_corr.for"
           enddo
           enddo
           enddo

        endif


      ENDIF

      end
