               l1 = l+1
               l2 = l-1
               l3 = l+incj
               l4 = l-incj
               l5 = l+inck
               l6 = l-inck

               xvol = 0.5/vol(lvo)

               u1 = (rop(l,2)+rop(l2,2))
               u2 = (rop(l,2)+rop(l3,2))
               u3 = (rop(l,2)+rop(l1,2))
               u4 = (rop(l,2)+rop(l4,2))
               u5 = (rop(l,2)+rop(l5,2))
               u6 = (rop(l,2)+rop(l6,2))
               !dudx
               dudx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx  ) 

               auijuij = dudx*dudx

               !dudy
               rotz =-(   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 
               dudy = rotz
               auijuij = auijuij + rotz*rotz

               !dudz
               roty = (u5 - u6)*tkz
               dudz = roty
               auijuij = auijuij + roty*roty


               u1 = (rop(l,3)+rop(l2,3))
               u2 = (rop(l,3)+rop(l3,3))
               u3 = (rop(l,3)+rop(l1,3))
               u4 = (rop(l,3)+rop(l4,3))
               u5 = (rop(l,3)+rop(l5,3))
               u6 = (rop(l,3)+rop(l6,3))
               !dvdx
               dvdx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx ) 

               rotz    = (rotz + dvdx)*xvol
               auijuij = auijuij + dvdx*dvdx

               !dvdy
               dvdy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 

               auijuij = auijuij + dvdy*dvdy

               !dvdzadding a source term to the momentum equation
               dvdz = (u5 - u6)*tkz

               rotx    =-dvdz 
               auijuij =  auijuij + dvdz*dvdz

               u1 = (rop(l,4)+rop(l2,4))
               u2 = (rop(l,4)+rop(l3,4))
               u3 = (rop(l,4)+rop(l1,4))
               u4 = (rop(l,4)+rop(l4,4))
               u5 = (rop(l,4)+rop(l5,4))
               u6 = (rop(l,4)+rop(l6,4))
               !dwdx
               dwdx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx )

               roty    =  (roty - dwdx)*xvol
               auijuij =  auijuij + dwdx*dwdx
               !dwdy
               dwdy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy  ) 

               rotx    = (rotx + dwdy)*xvol
               auijuij = auijuij + dwdy*dwdy

               !dwdz
               dwdz = (u5 - u6)*tkz

               auijuij =  auijuij + dwdz*dwdz

               !! mise a jour  auijuij par le volume
               rot     = sqrt(rotx*rotx+roty*roty+rotz*rotz)

               S11   = dudx
               S22   = dvdy
               S33   = dwdz
               S12   = 0.5*(dudy + dvdx)
               S13   = 0.5*(dudz + dwdx)
               S23   = 0.5*(dvdz + dwdy)
               St    = S11**2+S22**2+S33**2+2*S12**2+2*S13**2+2*S23**2
               St    = sqrt(2.*St) * xvol

               auijuij   = auijuij*xvol*xvol
               auijuij   = max( auijuij, 1.e-27)

               !formulation compressible complete
               u1 = (rop(l,6)+rop(l2,6))
               u2 = (rop(l,6)+rop(l3,6))
               u3 = (rop(l,6)+rop(l1,6))
               u4 = (rop(l,6)+rop(l4,6))
               u5 = (rop(l,6)+rop(l5,6))
               u6 = (rop(l,6)+rop(l6,6))
               !dudx
               dudx  = (   u3*tix1 - u1*tix
     &                   + u2*tjx1 - u4*tjx ) 

               !dudy
               dudy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy ) 
               !dudz
               dudz = (u5 - u6)*tkz


              u1 = rop(l,6)*rop(l,1) + rop(l2,6)*rop(l2,1)
              u2 = rop(l,6)*rop(l,1) + rop(l3,6)*rop(l3,1)
              u3 = rop(l,6)*rop(l,1) + rop(l1,6)*rop(l1,1)
              u4 = rop(l,6)*rop(l,1) + rop(l4,6)*rop(l4,1)
              u5 = rop(l,6)*rop(l,1) + rop(l5,6)*rop(l5,1)
              u6 = rop(l,6)*rop(l,1) + rop(l6,6)*rop(l6,1)
              !dudx
               dudx = dudx* (  u3*tix1 - u1*tix
     &                       + u2*tjx1 - u4*tjx )

              !dudy
               dudy = dudy* (  u3*tiy1 - u1*tiy
     &                       + u2*tjy1 - u4*tjy ) 

              !dudz
               dudz = dudz* (u5 - u6)*tkz

               anvisc    = c1*(dudx+dudy+dudz)*xvol*xvol

