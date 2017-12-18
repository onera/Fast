 
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
 
       b11=coe(l1,2)*signe
       b12=tcx
       b21=(tcx*ph2-u*qn)
       b22=(qn-tcx*gam2*u)+coe(l1,2)*signe
       b23=(-tcx*gamm1*v)
       b24=(-tcx*gamm1*w)
       b25=tcx*gamm1
       b31=(-v*qn)
       b32=(tcx*v)
       b33=(qn)+coe(l1,2)*signe
       b41=(-w*qn)
       b42=(tcx*w)
       b44=(qn)+coe(l1,2)*signe
       b51=qn*(ph2 - h)
       b52=(tcx*h-gamm1*u*qn)
       b53=(     -gamm1*v*qn)
       b54=(-gamm1*w*qn)
       b55=(gam1*qn)+coe(l1,2)*signe

      b6 = (qn+coe(l1,2)*signe)*drodm(l1,6)

      b1= b11*drodm(l1,1)+b12*drodm(l1,2)
      b2= b21*drodm(l1,1)+b22*drodm(l1,2)+b23*drodm(l1,3)
     1   +b24*drodm(l1,4)+b25*drodm(l1,5)
      b3= b31*drodm(l1,1)+b32*drodm(l1,2)+b33*drodm(l1,3)
      b4= b41*drodm(l1,1)+b42*drodm(l1,2)
     1   +b44*drodm(l1,4)
      b5= b51*drodm(l1,1)+b52*drodm(l1,2)+b53*drodm(l1,3)
     1   +b54*drodm(l1,4)+b55*drodm(l1,5)

      drodm(l,1) = drodm(l,1)+b1*xal 
      drodm(l,2) = drodm(l,2)+b2*xal 
      drodm(l,3) = drodm(l,3)+b3*xal 
      drodm(l,4) = drodm(l,4)+b4*xal 
      drodm(l,5) = drodm(l,5)+b5*xal 
      drodm(l,6) = drodm(l,6)+b6*xal
