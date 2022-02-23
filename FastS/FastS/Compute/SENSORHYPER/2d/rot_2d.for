               l1 = l+1
               l2 = l-1
               l3 = l+incj
               l4 = l-incj

               xvol = 0.5/vol(lvo)

               u1 = (rop(l,2)+rop(l2,2))
               u2 = (rop(l,2)+rop(l3,2))
               u3 = (rop(l,2)+rop(l1,2))
               u4 = (rop(l,2)+rop(l4,2))
               !dudx
               dudx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx  ) 

               divx = dudx*xvol
               !dudy
               rotz =-(   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 



               u1 = (rop(l,3)+rop(l2,3))
               u2 = (rop(l,3)+rop(l3,3))
               u3 = (rop(l,3)+rop(l1,3))
               u4 = (rop(l,3)+rop(l4,3))
               !dvdx
               dudx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx ) 

               rotz    = (rotz + dudx)*xvol

               !dvdy
               dudx = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 
               divy = dudx*xvol


               !! mise a jour  auijuij par le volume
               rot     = rotz*rotz
               div     = (divx + divy)*(divx + divy)
               umag    = sqrt(rop(l,2)**2+rop(l,3)**2+rop(l,4)**2.)
               mach    = min(umag/sqrt(gam*rgp*rop(l,5)),1.)**4
               voldes  = max(vol(lvo) , 1.e-27) 
               wig(l+wigd)=mach*(0.5+0.5*tanh((div /
     &          (div + (coef1*umag/(voldes**(1./3.)))**2. + 1d-14)
     &         *((div/(div+rot+1d-14)))-coef2)*40.*3.1415926535))
     
     
     
               !wig(l+wigd)=1.







