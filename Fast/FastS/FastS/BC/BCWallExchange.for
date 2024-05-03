            m     = lsample

            tcx  = tijk(lmtr,ic)*ci_mtr
            tcy  = tijk(lmtr,jc)*cj_mtr
            tcz  = tijk(lmtr,kc)*ck_mtr

            surf = sqrt(tcx*tcx+tcy*tcy+tcz*tcz)
            surf = max(surf,1e-30)

            nx = tcx/surf
            ny = tcy/surf
            nz = tcz/surf

            ventx      = ventijk(ldp ,1      )*c_ale
            venty      = ventijk(ldp ,2      )*c_ale
            ventz      = ventijk(ldp ,kc_vent)*ck_vent*c_ale

          qen =(  ventijk(ldp ,1      )*tcx
     &           +ventijk(ldp ,2      )*tcy
     &           +ventijk(ldp ,kc_vent)*tcz*ck_vent
     &         )*c_ale

            ! Calculate tangential velocity and its norm
            u = rop(m,2) 
            v = rop(m,3) 
            w = rop(m,4) 

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
            dist  = xmut(m, 2) !wall distance

            ! Wall function for the velocity proposed in
            ! "Explicit expression for the smooth wall velocity
            ! distribution in a turbulent boundary layer"
            ! AIAA Journal. Musker (1979)

            ! Initial prediction based on linear profile (u+=y+)
            if(sample.eq.1) then
              utau = sqrt( (xmut(ldjr,1)*unorm)/(rop(m,1)*dist) ) 
            else
              utau=((unorm**7)*xmut(ldjr,1)/(rop(m,1)*dist)/(8.3**7)
     &             )**0.125
            endif
c              !utau = sqrt( (xmut(ldjr,1)*unorm)/(rop(m,1)*dist) ) 

            ! Newton method to determine utau using Musker's wall law
            iter = 0 
            aa    = rop(m,1)*dist/xmut(ldjr,1)
            yplus = aa*utau
            g1    = (yplus+10.6)**9.6
            g2    = yplus*(yplus-8.15) + 86
            g3    = (2*yplus-8.15)/16.7
            f = 5.424*atan(g3) + log10(g1/(g2*g2)) - unorm/utau - 3.52

c           if(k.eq.6.and.i.eq.14) then
c           write(*,'(a,6f20.12,i3)') 'init',utau, yplus,unorm,xmut(l,1),
c     &           rop(m,1),dist,sample
c           endif

            do while (abs(f).gt.1e-6.and.iter.le.50)
 
             iter = iter + 1
 
             aa    = rop(m,1)*dist/xmut(ldjr,1)
             yplus = aa*utau
             g1    = (yplus+10.6)**9.6
             g2    = yplus*(yplus-8.15) + 86
             g3    = (2*yplus-8.15)/16.7
 
             tp =aa*(9.6*((yplus + 10.6)**8.6)- g1/g2*( 4*yplus-16.30) )
 
             fp = 0.649580838323*aa/(1 + g3*g3 )  +  tp/(g1*2.302585093)
     &        + unorm/(utau*utau)
 
             !blinder fp
             utau = abs(utau - f/fp)

             aa    = rop(m,1)*dist/xmut(ldjr,1)
             yplus = aa*utau
             g1    = (yplus+10.6)**9.6
             g2    = yplus*(yplus-8.15) + 86
             g3    = (2*yplus-8.15)/16.7

             f = 5.424*atan(g3) + log10(g1/(g2*g2)) - unorm/utau - 3.52
            end do

           ! test affinage utau en passant par psiroe
           !utau = utau*param_real(PSIROE)

           yplus = rop(ldjr,1)*utau*xmut(l,2)/xmut(ldjr,1)
           umod = utau*(5.424*atan((2.*yplus - 8.15)/16.7)
     &            + log10((yplus + 10.6)**9.6
     &            /(yplus**2 - 8.15*yplus + 86.)**2) - 3.52)

           !umod = abs(umod)*param_real(PSIROE)

           rop(l,1) = rop(m,1)
           rop(l,2) = umod*ut/unorm + xmut(l,2)/xmut(m,2)*un
           rop(l,3) = umod*vt/unorm + xmut(l,2)/xmut(m,2)*vn
           rop(l,4) = umod*wt/unorm + xmut(l,2)/xmut(m,2)*wn
           rop(l,5) = rop(m,5)

           pwall = rop(l,5)*rop(l,1)
           twall = rop(l,5) 
     &   + 0.5*prandlt**0.33333/(cv*gam)*(unorm**2-umod**2) ! Crocco-Busemann
           rop(l,5) = twall
           rop(l,1) = pwall/twall
