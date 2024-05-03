      !dir ijk
      !sciacovelli
      lw  = l + w_tg
      mach=sqrt(gam*rgp*rop(l,5))
      mach=10.0*(1.0/tc_tg)*(1.0/mach)*div
      mach=0.5*(1.0-tanh(2.5d0+mach))
              
      !jameson
      p2 = rop(l-inc_tg,1)*rop(l-inc_tg,5)
      p0 = rop(l       ,1)*rop(l       ,5)
      p1 = rop(l+inc_tg,1)*rop(l+inc_tg,5)
      wig(lw)= abs( p2 - 2.*p0 + p1)/abs( p2 + 2.*p0 + p1)
      !assemblage
      !div=max(abs(div),0.00001)
      !rot=max(0.00001,rot)
      !mach = 0.1
      wig(lw)= min( 2.0*div*div/(div*div+rot+1d-15)*wig(lw)*mach, 1.0)
      !wig(lw)= max( wig(lw), 1.e-15)
      wig(lw)= 1.0 - wig(lw)

