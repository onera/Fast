 
       l1  = l -inci

c-----Metrque au centre en + ou - 1 : 

       r     = rop(l1,1)
       u     = rop(l1,2)
       v     = rop(l1,3)
       w     = rop(l1,4)
       t     = rop(l1,5)
       q2    = .5*(u*u+v*v+w*w)
       h     = cp*t + q2
       ph2   = gamm1*q2

       qn    = tcx*u
       qen   = tcx*venti(lv)
 
       b11=coe(l1,2)*signe -qen
       b12=tcx
       b21=(tcx*ph2-u*qn)
       b22=(qn-tcx*gam2*u)+coe(l1,2)*signe -qen
       b23=(-tcx*gamm1*v)
       b24=(-tcx*gamm1*w)
       b25=tcx*gamm1
       b31=(-v*qn)
       b32=(tcx*v)
       b33=(qn)+coe(l1,2)*signe -qen
       b41=(-w*qn)
       b42=(tcx*w)
       b44=(qn)+coe(l1,2)*signe -qen
       b51=qn*(ph2 - h)
       b52=(tcx*h-gamm1*u*qn)
       b53=(     -gamm1*v*qn)
       b54=(-gamm1*w*qn)
       b55=(gam1*qn)+coe(l1,2)*signe -qen

      b6 = (qn+coe(l1,2)*signe -qen)*drodm_out(l1,6)

      b1= b11*drodm_out(l1,1)+b12*drodm_out(l1,2)
      b2= b21*drodm_out(l1,1)+b22*drodm_out(l1,2)+b23*drodm_out(l1,3)
     1   +b24*drodm_out(l1,4)+b25*drodm_out(l1,5)
      b3= b31*drodm_out(l1,1)+b32*drodm_out(l1,2)+b33*drodm_out(l1,3)
      b4= b41*drodm_out(l1,1)+b42*drodm_out(l1,2)
     1   +b44*drodm_out(l1,4)
      b5= b51*drodm_out(l1,1)+b52*drodm_out(l1,2)+b53*drodm_out(l1,3)
     1   +b54*drodm_out(l1,4)+b55*drodm_out(l1,5)
