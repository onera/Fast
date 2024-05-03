      l1  = l -inck
      lv1 = lv-inck_ven
 
      tcz = tk(lt,1)

      r    = rop(l1,1)
      u    = rop(l1,2)
      v    = rop(l1,3)
      w    = rop(l1,4)
      t    = rop(l1,5)
      q2   = .5*(u*u+v*v+w*w)
      h    = cp*t + q2
      ph2  = gamm1*q2
      
      qn    = tcz*w
      qen = tcz*ventk(lv1+v3ven)*ck_vent

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
 
      b6 = (qn+coe(l1,4)*signe -qen)*drodm_out(l1,6)

      b1= b11*drodm_out(l1,1)                   + b14*drodm_out(l1,4)
      b2= b21*drodm_out(l1,1)+ b22*drodm_out(l1,2) + b24*drodm_out(l1,4)
      b3= b31*drodm_out(l1,1)+ b33*drodm_out(l1,3) + b34*drodm_out(l1,4)
      b4= b41*drodm_out(l1,1)+ b42*drodm_out(l1,2) + b43*drodm_out(l1,3)
     &   +b44*drodm_out(l1,4)+ b45*drodm_out(l1,5)
      b5= b51*drodm_out(l1,1)+ b52*drodm_out(l1,2) + b53*drodm_out(l1,3)
     &   +b54*drodm_out(l1,4)+ b55*drodm_out(l1,5)
