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
      dudx = ( u3 - u1)*tix  !dudx
      rotz =-( u2 - u4)*tjy  !dudy
      roty = ( u5 - u6)*tkz  !dudz
      div = dudx*xvol


      u1 = (rop(l,3)+rop(l2,3))
      u2 = (rop(l,3)+rop(l3,3))
      u3 = (rop(l,3)+rop(l1,3))
      u4 = (rop(l,3)+rop(l4,3))
      u5 = (rop(l,3)+rop(l5,3))
      u6 = (rop(l,3)+rop(l6,3))
      
      dudx = ( u3 - u1)*tix   !dvdx
      rotz  = (rotz + dudx)*xvol

      dudx = ( u2 - u4)*tjy   !dvdy
      div  = div + dudx*xvol

      dudx = ( u5 - u6)*tkz !dvdz
      rotx =-dudx 

      u1 = (rop(l,4)+rop(l2,4))
      u2 = (rop(l,4)+rop(l3,4))
      u3 = (rop(l,4)+rop(l1,4))
      u4 = (rop(l,4)+rop(l4,4))
      u5 = (rop(l,4)+rop(l5,4))
      u6 = (rop(l,4)+rop(l6,4))
      
      dudx = ( u3 - u1)*tix  !dwdx
      roty = (roty - dudx)*xvol
      dudx = ( u2 - u4)*tjy  !dwdy
      rotx = (rotx + dudx)*xvol
      dudx = ( u5 - u6)*tkz  !dwdz
      div= div + dudx*xvol


      !! mise a jour  auijuij par le volume
      rot     = rotx*rotx+roty*roty+rotz*rotz
