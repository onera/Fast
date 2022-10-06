               !calcul de delta et gradient vitesse
               tix  = ti(lt,1)
               tiy  = ti(lt,2)
               tix1 = ti(lt+inci_mtr,1)
               tiy1 = ti(lt+inci_mtr,2)

               tjx  = tj(lt,1)
               tjy  = tj(lt,2)
               tjx1 = tj(lt+incj_mtr,1)
               tjy1 = tj(lt+incj_mtr,2)
c
               l1 = l+1
               l2 = l-1
               l3 = l+incj
               l4 = l-incj

               u1 = (rop(l,2)+rop(l2,2))
               u2 = (rop(l,2)+rop(l3,2))
               u3 = (rop(l,2)+rop(l1,2))
               u4 = (rop(l,2)+rop(l4,2))

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

               rotz    = rotz + dudx

               !! mise a jour rot et auijuij par le volume
               xvol = 0.5/vol(lvo)
               rot     = abs(rotz)* xvol

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
