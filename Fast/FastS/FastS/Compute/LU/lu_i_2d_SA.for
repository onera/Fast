 
       l1  = l -inci
       l1s = ls-incis
       lt1 = lt-inci_mtr
       lt2 = lt-inci2_mtr
 
       tcx = 0.5*(ti(lt2,1)+ti(lt1,1))
       tcy = 0.5*(ti(lt2,2)+ti(lt1,2))

       r    = rop(l1,1)
       u    = rop(l1,2)
       v    = rop(l1,3)
       t    = rop(l1,5)
       q2   = .5*(u*u+v*v)
       h    = cp*t + q2
       ph2  = gamm1*q2

       qn    = tcx*u+tcy*v
 
       b11=coe(l1,2)*signe
       b12=tcx
       b13=tcy
       b21=(tcx*ph2-u*qn)
       b22=(qn-tcx*gam2*u)+coe(l1,2)*signe
       b23=(tcy*u-tcx*gamm1*v)
       b25=tcx*gamm1
       b31=(tcy*ph2-v*qn)
       b32=(tcx*v-tcy*gamm1*u)
       b33=(qn-tcy*gam2*v)+coe(l1,2)*signe
       b35=tcy*gamm1
       b51=qn*(ph2 - h)
       b52=(tcx*h-gamm1*u*qn)
       b53=(tcy*h-gamm1*v*qn)
       b55=(gam1*qn)+coe(l1,2)*signe

      b6 = (qn+coe(l1,2)*signe)*drodm_out(l1s,6)*1.e0

      b1= b11*drodm_out(l1s,1)+b12*drodm_out(l1s,2)+b13*drodm_out(l1s,3)
      b2= b21*drodm_out(l1s,1)+b22*drodm_out(l1s,2)+b23*drodm_out(l1s,3)
     1   +b25*drodm_out(l1s,5)
      b3= b31*drodm_out(l1s,1)+b32*drodm_out(l1s,2)+b33*drodm_out(l1s,3)
     1   +b35*drodm_out(l1s,5)
      b5= b51*drodm_out(l1s,1)+b52*drodm_out(l1s,2)+b53*drodm_out(l1s,3)
     1   +b55*drodm_out(l1s,5)


