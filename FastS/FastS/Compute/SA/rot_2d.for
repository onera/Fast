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

               auijuij = dudx*dudx

               !dudy
               rotz =-(   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 
               dudy = rotz
               auijuij = auijuij + rotz*rotz


               u1 = (rop(l,3)+rop(l2,3))
               u2 = (rop(l,3)+rop(l3,3))
               u3 = (rop(l,3)+rop(l1,3))
               u4 = (rop(l,3)+rop(l4,3))
               !dvdx
               dvdx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx ) 

               rotz    = (rotz + dvdx)*xvol
               auijuij = auijuij + dvdx*dvdx

               !dvdy
               dvdy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 

               auijuij = auijuij + dvdy*dvdy


               !! mise a jour  auijuij par le volume
               rot     = abs(rotz)

               S11   = dudx
               S22   = dvdy
               S12   = 0.5*(dudy + dvdx)
               St    = S11**2+S22**2+2*S12**2
               St    = sqrt(2*St) * xvol


               auijuij   = auijuij*xvol*xvol
               auijuij   = max( auijuij, 1.e-27)

               !formulation compressible complete
               u1 = (rop(l,6)+rop(l2,6))
               u2 = (rop(l,6)+rop(l3,6))
               u3 = (rop(l,6)+rop(l1,6))
               u4 = (rop(l,6)+rop(l4,6))
               !dudx
               dudx  = (   u3*tix1 - u1*tix
     &                   + u2*tjx1 - u4*tjx ) 

               !dudy
               dudy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 


              u1 = rop(l,6)*rop(l,1) + rop(l2,6)*rop(l2,1)
              u2 = rop(l,6)*rop(l,1) + rop(l3,6)*rop(l3,1)
              u3 = rop(l,6)*rop(l,1) + rop(l1,6)*rop(l1,1)
              u4 = rop(l,6)*rop(l,1) + rop(l4,6)*rop(l4,1)
              !dudx
               dudx = dudx* (  u3*tix1 - u1*tix
     &                       + u2*tjx1 - u4*tjx )

              !dudy
               dudy = dudy* (  u3*tiy1 - u1*tiy
     &                       + u2*tjy1 - u4*tjy ) 

               anvisc    = c1*(dudx+dudy)*xvol*xvol

