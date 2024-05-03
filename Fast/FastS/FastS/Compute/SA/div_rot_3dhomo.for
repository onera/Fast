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
      dudx = (   u3*tix1 - u1*tix + u2*tjx1 - u4*tjx  ) 
      div = dudx*xvol

      !dudy
      rotz =-(   u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy ) 

      !dudz
      roty = (u5 - u6)*tkz

      u1 = (rop(l,3)+rop(l2,3))
      u2 = (rop(l,3)+rop(l3,3))
      u3 = (rop(l,3)+rop(l1,3))
      u4 = (rop(l,3)+rop(l4,3))
      u5 = (rop(l,3)+rop(l5,3))
      u6 = (rop(l,3)+rop(l6,3))
      !dvdx
      dudx = (   u3*tix1 - u1*tix + u2*tjx1 - u4*tjx ) 

      rotz    = (rotz + dudx)*xvol

      !dvdy
      dudx = (   u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy ) 

      div = div + dudx*xvol

      !dvdzadding a source term to the momentum equation
      dudx = (u5 - u6)*tkz

      rotx    =-dudx 

      u1 = (rop(l,4)+rop(l2,4))
      u2 = (rop(l,4)+rop(l3,4))
      u3 = (rop(l,4)+rop(l1,4))
      u4 = (rop(l,4)+rop(l4,4))
      u5 = (rop(l,4)+rop(l5,4))
      u6 = (rop(l,4)+rop(l6,4))
      !dwdx
      dudx = (   u3*tix1 - u1*tix + u2*tjx1 - u4*tjx )

      roty    =  (roty - dudx)*xvol
      !dwdy
      dudx = (   u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy ) 

      rotx    = (rotx + dudx)*xvol

      !dwdz
      dudx = (u5 - u6)*tkz
      div= div + dudx*xvol


      !! mise a jour  auijuij par le volume
      rot     = rotx*rotx+roty*roty+rotz*rotz
