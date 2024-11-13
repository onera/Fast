C... End of Newton and building of the variables
          vng  = (wng+qen)
C...      Absolute velocity module
          vg   = vng*usd0n
          
          b    = 1. - (vng**2)*usd0n2/(2.*ha(li))
          b = max(0.2,b) !cutoff robustesse

          pg   = pa(li)*b**gam2
          rog  = gam2*pg/(ha(li)*b)
          
          Tg   = pg/(rog*rgp)

          rop(l,1) = rog
          rop(l,2) = vg*d0x(li)
          rop(l,3) = vg*d0y(li)
          rop(l,4) = vg*d0z(li)
          rop(l,5) = Tg
