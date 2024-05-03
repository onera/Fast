      !dir i
      !sciacovelli
      lw  = l + wigi
      mach=sqrt(gam*rgp*rop(l,5))
      mach=10.0*(1.0/tci)*(1.0/mach)*div
      mach=0.5*(1.0-tanh(2.5d0+mach))
              
      !jameson
      p2 = rop(l2,1)*rop(l2,5)
      p0 = rop(l ,1)*rop(l ,5)
      p1 = rop(l1,1)*rop(l1,5)
      wig(lw)= abs( p2 - 2.*p0 + p1)/abs( p2 + 2.*p0 + p1)
      !assemblage
      wig(lw)= min( 2.0*div*div/(div*div+rot+1d-15)*wig(lw)*mach, 1.0)
      wig(lw)= 1.0 - wig(lw)

      !dir j
      !sciacovelli
      lw  = l + wigj
      mach=sqrt(gam*rgp*rop(l,5))
      mach=10.0*(1.0/tcj)*(1.0/mach)*div
      mach=0.5*(1.0-tanh(2.5d0+mach))
      !jameson
      p4 = rop(l4,1)*rop(l4,5)
      p3 = rop(l3,1)*rop(l3,5)
      wig(lw)= abs( p4 - 2.*p0 + p3)/abs( p4 + 2.*p0 + p3)
      !assemblage
      wig(lw)= min( 2.0*div*div/(div*div+rot+1d-15)*wig(lw)*mach, 1.0)
      wig(lw)= 1.0 - wig(lw)
