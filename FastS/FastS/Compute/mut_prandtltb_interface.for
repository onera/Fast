      xmulam = coesut*sqrt(rop(il +v5)*temp01) / (1.+cmus1/rop(il +v5))

      xmulam = xmulam + coesut *sqrt( rop(ir +v5)*temp01) 
     &                       /(1.+cmus1/rop(ir +v5))

      xmutot = xmut(ir)+xmut(il)
      xmutur = xmutot - xmulam
      xktvol = ( xmulam*gam3 + xmutur*gam4)*volinv

      xmutvol=xmutot*volinv
