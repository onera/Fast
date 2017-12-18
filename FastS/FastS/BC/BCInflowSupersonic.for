      rop(l,1) =  ro*c3 + rop(l1,1)*c1  +rop(l2,1)*c2 

      !if(rop(l,1).le.0.3*ro) then  ! protection pour eviter valeur negative de ro et t
      if( abs(rop(l,1)/ro -1.).ge.0.7) then  ! protection pour eviter valeur negative de ro et t
        rop(l,1) = ro
        rop(l,2) = u
        rop(l,3) = v
        rop(l,4) = w
        rop(l,5) = t
      else
        rop(l,2) =   u*c3 + rop(l1,2)*c1  +rop(l2,2)*c2 
        rop(l,3) =   v*c3 + rop(l1,3)*c1  +rop(l2,3)*c2 
        rop(l,4) =   w*c3 + rop(l1,4)*c1  +rop(l2,4)*c2 
        rop(l,5) =   t*c3 + rop(l1,5)*c1  +rop(l2,5)*c2 
      endif
