               l1 = l+1
               l2 = l-1
               l3 = l+incj
               l4 = l-incj
               l5 = l+inck
               l6 = l-inck

               xvol = 0.5/vol(lvo)

               !!!VelocityX
               u1 = (rop(l +v2)+rop(l2 +v2))
               u2 = (rop(l +v2)+rop(l3 +v2))
               u3 = (rop(l +v2)+rop(l1 +v2))
               u4 = (rop(l +v2)+rop(l4 +v2))
               u5 = (rop(l +v2)+rop(l5 +v2))
               u6 = (rop(l +v2)+rop(l6 +v2))

               dudx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx 
     &                  + u5*tkx1 - u6*tkx  ) 

               dudy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy 
     &                  + u5*tky1 - u6*tky )

               dudz = (   u3*tiz1 - u1*tiz
     &                  + u2*tjz1 - u4*tjz 
     &                  + u5*tkz1 - u6*tkz ) 
               !!!VelocityY
               u1 = (rop(l +v3)+rop(l2 +v3))
               u2 = (rop(l +v3)+rop(l3 +v3))
               u3 = (rop(l +v3)+rop(l1 +v3))
               u4 = (rop(l +v3)+rop(l4 +v3))
               u5 = (rop(l +v3)+rop(l5 +v3))
               u6 = (rop(l +v3)+rop(l6 +v3))

               dvdx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx 
     &                  + u5*tkx1 - u6*tkx  ) 

               dvdy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy 
     &                  + u5*tky1 - u6*tky  ) 

               dvdz = (   u3*tiz1 - u1*tiz
     &                  + u2*tjz1 - u4*tjz 
     &                  + u5*tkz1 - u6*tkz  ) 
               !!!VelocityZ
               u1 = (rop(l +v4)+rop(l2 +v4))
               u2 = (rop(l +v4)+rop(l3 +v4))
               u3 = (rop(l +v4)+rop(l1 +v4))
               u4 = (rop(l +v4)+rop(l4 +v4))
               u5 = (rop(l +v4)+rop(l5 +v4))
               u6 = (rop(l +v4)+rop(l6 +v4))

               dwdx = (   u3*tix1 - u1*tix
     &                  + u2*tjx1 - u4*tjx 
     &                  + u5*tkx1 - u6*tkx  )

               dwdy = (   u3*tiy1 - u1*tiy
     &                  + u2*tjy1 - u4*tjy 
     &                  + u5*tky1 - u6*tky  ) 


               dwdz = (   u3*tiz1 - u1*tiz
     &                  + u2*tjz1 - u4*tjz 
     &                  + u5*tkz1 - u6*tkz  ) 

               rot( l + v1) = ( dwdy - dvdz )*xvol
               rot( l + v2) = ( dudz - dwdx )*xvol
               rot( l + v3) = ( dvdx - dudy )*xvol

               xmut(l) = (dudx*dudx + dvdy*dvdy + dwdz*dwdz)*2.
               xmut(l) = xmut(l) + (dvdx+dudy)*(dvdx+dudy)
               xmut(l) = xmut(l) + (dwdx+dudz)*(dwdx+dudz)
               xmut(l) = xmut(l) + (dwdy+dvdz)*(dwdy+dvdz)
               xmut(l) = xmut(l)*0.25                      !0.5"vol"*0.5"vol" 
