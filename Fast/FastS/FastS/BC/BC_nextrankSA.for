      rop(l,1) =  rop(l0,1)*c3 + rop(l1,1)*c1  +rop(l2,1)*c2 

      !if( abs(rop(l,1)/rop(l0,1) -1.).ge.0.7) then  

      if(rop(l,1).le.0.3*rop(l0,1)) then  ! protection pour eviter valeur negative de ro et t
        rop(l,1) = rop(l0,1)
        rop(l,2) = rop(l0,2)
        rop(l,3) = rop(l0,3)
        rop(l,4) = rop(l0,4)
        rop(l,5) = rop(l0,5)
        rop(l,6) = rop(l0,6)
      else
        rop(l,2) =   rop(l0,2)*c3 + rop(l1,2)*c1  +rop(l2,2)*c2 
        rop(l,3) =   rop(l0,3)*c3 + rop(l1,3)*c1  +rop(l2,3)*c2 
        rop(l,4) =   rop(l0,4)*c3 + rop(l1,4)*c1  +rop(l2,4)*c2 
        rop(l,5) =   rop(l0,5)*c3 + rop(l1,5)*c1  +rop(l2,5)*c2 
        rop(l,6) =   rop(l0,6)*c3 + rop(l1,6)*c1  +rop(l2,6)*c2 
      endif
