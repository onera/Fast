      l1  = l -inck
 
      r    = rop(l1,1)
      u    = rop(l1,2)
      v    = rop(l1,3)
      w    = rop(l1,4)
      t    = rop(l1,5)
      q2   = .5*(u*u+v*v+w*w)
      h    = cp*t + q2
      ph2  = gamm1*q2
      
      qn    = tcz*w
      qen   = tcz*ventk(lv+v3ven)*ck_vent

      b11= coe(l1,4)*signe -qen
      b14= tcz
      b21= -u*qn
      b22= qn+coe(l1,4)*signe -qen
      b24= tcz*u
      b31= -v*qn
      b33= qn+coe(l1,4)*signe -qen
      b34= tcz*v
      b41= tcz*ph2-w*qn
      b42=-tcz*gamm1*u
      b43=-tcz*gamm1*v
      b44= (qn-tcz*gam2*w)+coe(l1,4)*signe -qen
      b45= tcz*gamm1
      b51= qn*(ph2 - h)
      b52=-gamm1*u*qn
      b53=-gamm1*v*qn
      b54= tcz*h-gamm1*w*qn
      b55= gam1*qn+coe(l1,4)*signe -qen
 
      b6 = (qn+coe(l1,4)*signe -qen)*drodm(l1,6)

      b1= b11*drodm(l1,1)                   + b14*drodm(l1,4)
      b2= b21*drodm(l1,1) + b22*drodm(l1,2) + b24*drodm(l1,4)
      b3= b31*drodm(l1,1) + b33*drodm(l1,3) + b34*drodm(l1,4)
      b4= b41*drodm(l1,1) + b42*drodm(l1,2) + b43*drodm(l1,3)
     &   +b44*drodm(l1,4) + b45*drodm(l1,5)
      b5= b51*drodm(l1,1) + b52*drodm(l1,2) + b53*drodm(l1,3)
     &   +b54*drodm(l1,4) + b55*drodm(l1,5)

      drodm(l,1) = drodm(l,1)+b1*xal 
      drodm(l,2) = drodm(l,2)+b2*xal 
      drodm(l,3) = drodm(l,3)+b3*xal 
      drodm(l,4) = drodm(l,4)+b4*xal 
      drodm(l,5) = drodm(l,5)+b5*xal
      drodm(l,6) = drodm(l,6)+b6*xal
