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
      dudx = (   u3*tix1 - u1*tix + u2*tjx1 - u4*tjx  ) 
      div = dudx*xvol

      !dudy
      rotz =-(   u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy ) 

      u1 = (rop(l,3)+rop(l2,3))
      u2 = (rop(l,3)+rop(l3,3))
      u3 = (rop(l,3)+rop(l1,3))
      u4 = (rop(l,3)+rop(l4,3))
      !dvdx
      dudx = (   u3*tix1 - u1*tix + u2*tjx1 - u4*tjx ) 

      rotz    = (rotz + dudx)*xvol

      !dvdy
      dudx = (   u3*tiy1 - u1*tiy + u2*tjy1 - u4*tjy ) 
      div = div + dudx*xvol

      !! mise a jour  auijuij par le volume
      rot     = rotz*rotz
