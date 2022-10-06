               !calcul de delta et gradient vitesse
               tix  = ti(lt,1)
               tiy  = ti(lt,2)
               tix1 = ti(lt+inci_mtr,1)
               tiy1 = ti(lt+inci_mtr,2)

               tjx  = tj(lt,1)
               tjy  = tj(lt,2)
               tjx1 = tj(lt+incj_mtr,1)
               tjy1 = tj(lt+incj_mtr,2)

               tkz  = tk(lt,1)

               l1 = l+1
               l2 = l-1
               l3 = l+incj
               l4 = l-incj
               l5 = l+inck
               l6 = l-inck

               u1 = (rop(l,2)+rop(l2,2))
               u2 = (rop(l,2)+rop(l3,2))
               u3 = (rop(l,2)+rop(l1,2))
               u4 = (rop(l,2)+rop(l4,2))
               u5 = (rop(l,2)+rop(l5,2))
               u6 = (rop(l,2)+rop(l6,2))

               !dudx
               dudx = ( u3*tix1 - u1*tix + u2*tjx1 - u4*tjx )

               !dudy
               dudy = ( u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy )

               !dudz
               dudz = ( u5 - u6)*tkz 

               u1 = (rop(l,3)+rop(l2,3))
               u2 = (rop(l,3)+rop(l3,3))
               u3 = (rop(l,3)+rop(l1,3))
               u4 = (rop(l,3)+rop(l4,3))
               u5 = (rop(l,3)+rop(l5,3))
               u6 = (rop(l,3)+rop(l6,3))

               !dvdx
               dvdx = (  u3*tix1 - u1*tix + u2*tjx1 - u4*tjx )

               !dvdy
               dvdy = (  u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy )

               !dvdz
               dvdz = (  u5 - u6)*tkz 

               u1 = (rop(l,4)+rop(l2,4))
               u2 = (rop(l,4)+rop(l3,4))
               u3 = (rop(l,4)+rop(l1,4))
               u4 = (rop(l,4)+rop(l4,4))
               u5 = (rop(l,4)+rop(l5,4))
               u6 = (rop(l,4)+rop(l6,4))

               !dwdx
               dwdx = (u3*tix1 - u1*tix + u2*tjx1 - u4*tjx )

               !dwdy
               dwdy = (u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy )

               !dwdz
               dwdz = (u5 - u6)*tkz 

               rotx    = dwdy - dvdz
               roty    = dudz - dwdx
               rotz    = dvdx - dudy

               !! mise a jour rot et auijuij par le volume
               xvol = 0.5/vol(lvo)
               rot     = sqrt(rotx*rotx+roty*roty+rotz*rotz)* xvol

               S11   = dudx
               S22   = dvdy
               S33   = dwdz
               S12   = 0.5*(dudy + dvdx)
               S13   = 0.5*(dudz + dwdx)
               S23   = 0.5*(dvdz + dwdy)
               St    = S11**2+S22**2+S33**2+2*S12**2+2*S13**2+2*S23**2
               St    = sqrt(2*St) * xvol

               cc = -3.5*rop(l,1)*rop(l,6)**2/(gam*rgp*rop(l,5))
     &              *(dudx**2+dudy**2+dudz**2
     &               +dvdx**2+dvdy**2+dvdz**2
     &               +dwdx**2+dwdy**2+dwdz**2)

               !formulation compressible complete
               u1 = (rop(l,6)+rop(l2,6))
               u2 = (rop(l,6)+rop(l3,6))
               u3 = (rop(l,6)+rop(l1,6))
               u4 = (rop(l,6)+rop(l4,6))
               u5 = (rop(l,6)+rop(l5,6))
               u6 = (rop(l,6)+rop(l6,6))
               !dudx
               dudx = (   u3*tix1 - u1*tix + u2*tjx1 - u4*tjx )
               !dudy
               dudy = (   u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy )
               !dudz
               dudz = (   u5 - u6)*tkz 

              u1 = rop(l,6)*rop(l,1) + rop(l2,6)*rop(l2,1)
              u2 = rop(l,6)*rop(l,1) + rop(l3,6)*rop(l3,1)
              u3 = rop(l,6)*rop(l,1) + rop(l1,6)*rop(l1,1)
              u4 = rop(l,6)*rop(l,1) + rop(l4,6)*rop(l4,1)
              u5 = rop(l,6)*rop(l,1) + rop(l5,6)*rop(l5,1)
              u6 = rop(l,6)*rop(l,1) + rop(l6,6)*rop(l6,1)
              !dudx
              dudx = dudx* (  u3*tix1 - u1*tix + u2*tjx1 - u4*tjx )

              !dudy
              dudy = dudy* (  u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy )

              !dudz
              dudz = dudz* (  u5 - u6)*tkz  

              anvisc    = ( c1*(dudx+dudy+dudz) + cc ) *xvol*xvol
