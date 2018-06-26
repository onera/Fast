       l1  = l -incj
c-----Metrque au centre en + ou - 1 : 

       r     = rop(l1,1)
       u     = rop(l1,2)
       v     = rop(l1,3)
       w     = rop(l1,4)
       t     = rop(l1,5)
       q2    = .5*(u*u+v*v+w*w)
       h     = cp*t + q2
       ph2   = gamm1*q2

       qn    = tcy*v
       qen   = tcy*ventj(lv+v2ven)
 
 
       b11=coe(l1,3)*signe -qen
       b13=tcy
       b21=(-u*qn)
       b22=qn+coe(l1,3)*signe -qen
       b23=tcy*u
       b31=(tcy*ph2-v*qn)
       b32=(-tcy*gamm1*u)
       b33=(qn-tcy*gam2*v)+coe(l1,3)*signe -qen
       b34=(-tcy*gamm1*w)
       b35=tcy*gamm1
       b41=(-w*qn)
       b43=(tcy*w)
       b44=(qn)+coe(l1,3)*signe -qen
       b51=qn*(ph2 - h)
       b52=(     -gamm1*u*qn)
       b53=(tcy*h-gamm1*v*qn)
       b54=(-gamm1*w*qn)
       b55=(gam1*qn)+coe(l1,3)*signe -qen
 

      b1= b11*drodm_out(l1,1)                +b13*drodm_out(l1,3)
      b2= b21*drodm_out(l1,1)+b22*drodm_out(l1,2)+b23*drodm_out(l1,3)
      b3= b31*drodm_out(l1,1)+b32*drodm_out(l1,2)+b33*drodm_out(l1,3)
     1   +b34*drodm_out(l1,4)+b35*drodm_out(l1,5)
      b4= b41*drodm_out(l1,1)                +b43*drodm_out(l1,3)
     1   +b44*drodm_out(l1,4)
      b5= b51*drodm_out(l1,1)+b52*drodm_out(l1,2)+b53*drodm_out(l1,3)
     1   +b54*drodm_out(l1,4)+b55*drodm_out(l1,5)
