            tcx = tijk(lt +vmtr*ic)*ci_mtr !surface normal x
            tcy = tijk(lt +vmtr*jc)*cj_mtr !surface normal y
            tcz = tijk(lt +vmtr*kc)*ck_mtr !surface normal z

            surf = sqrt(tcx*tcx+tcy*tcy+tcz*tcz) ! magnitude of surface normal
            surf = max(surf,1e-30)

            nx = tcx/surf ! normalized surface normal x
            ny = tcy/surf ! normalized surface normal y
            nz = tcz/surf ! normalized surface normal z

            ventx= ventijk(lven              )*c_ale*mobile_coef          ! rotational/translational velocity
            venty= ventijk(lven +vven        )*c_ale*mobile_coef          ! rotational/translational velocity 
            ventz= ventijk(lven +vven*kc_vent)*ck_vent*c_ale*mobile_coef  ! rotational/translational velocity       

          qen =(  ventijk(lven              )*tcx
     &           +ventijk(lven +vven        )*tcy
     &           +ventijk(lven +vven*kc_vent)*tcz*ck_vent
     &         )*c_ale     ! vent (dot) surface normal

            ! Calculate velocity repere tournant (eventuel) 
            u = rop(m +v2) -ventx
            v = rop(m +v3) -venty
            w = rop(m +v4) -ventz

            u_int = nx*u + ny*v + nz*w ! normal velocity au point sample

            un = u_int*nx
            vn = u_int*ny
            wn = u_int*nz

            ! Calculate tangential velocity and its norm
            ut = u-un
            vt = v-vn
            wt = w-wn
            unorm = sqrt(ut*ut+vt*vt+wt*wt)

            r =0.5*(rop(l+v1)+rop(iadrf+v1))
            p =0.5*(rop(l+v5)*rop(l+v1)+rop(iadrf+v5)*rop(iadrf+v1))*rgp

            dist = xmut(m+ v2) !wall distance

            ! Wall function for the velocity proposed in
            ! "Explicit expression for the smooth wall velocity
            ! distribution in a turbulent boundary layer"
            ! AIAA Journal. Musker (1979)

            !Visco wall
            t1    = 0.5*(rop(l+v5)+rop(iadrf+v5))
            mu_w  = coesut * sqrt(t1)*t1/(t1+cmus1)
            rho_w = 0.5*(rop(l+v1)+rop(iadrf+v1))

            ! Initial prediction based on linear profile (u+=y+)
            if(sampling.eq.1) then
              utau = sqrt( (mu_w*unorm)/(rho_w*dist) ) 
            else
              utau=((unorm**7)*mu_w/(rho_w*dist)/(8.3**7))**0.125
            endif

            aa    = rho_w*dist/mu_w

            ! Newton method to determine utau using Musker's wall law
            iter = 0 
            yplus = aa*utau
            l1    = (yplus+10.6)**9.6
            l2    = yplus*(yplus-8.15) + 86
            l3    = (2*yplus-8.15)/16.7
            f = 5.424*atan(l3) + log10(l1/(l2*l2)) - unorm/utau - 3.52

c           if(k.eq.1.and.i.eq.100) then
c       write(*,'(a,6f20.12,i3)') 'utau y+ unorm mu, rho, dist',
c     &    utau,yplus,unorm,xmut(ldjr+v1),
c     &           rop(m +v1),dist,sampling
c           endif

            do while (abs(f).gt.1e-6.and.iter.le.50)
 
             iter = iter + 1
 
             yplus = aa*utau
             l1    = (yplus+10.6)**9.6
             l2    = yplus*(yplus-8.15) + 86
             l3    = (2*yplus-8.15)/16.7
 
             tp =aa*(9.6*((yplus + 10.6)**8.6)- l1/l2*( 4*yplus-16.30) )
 
             fp = 0.649580838323*aa/(1 + l3*l3 )  +  tp/(l1*2.302585093)
     &        + unorm/(utau*utau)
 
             !blinder fp
             utau = abs(utau - f/fp)

             yplus = aa*utau
             l1    = (yplus+10.6)**9.6
             l2    = yplus*(yplus-8.15) + 86
             l3    = (2*yplus-8.15)/16.7

             f = 5.424*atan(l3) + log10(l1/(l2*l2)) - unorm/utau - 3.52
            end do

            utau2 = utau**2
            tauw = r*utau2

            u     = 0.5*(rop(l+v2)+rop(iadrf+v2))
            v     = 0.5*(rop(l+v3)+rop(iadrf+v3))
            w     = 0.5*(rop(l+v4)+rop(iadrf+v4))

            ! determination vitesse normale interface
            ! relative for rotational/translational
            !u_int= tcx*u +tcy*v +tcz*w -qen 
            u_int= 0.
