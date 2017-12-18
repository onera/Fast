CC    IN : rop
CC    OUT: roi,ui,vi,wi,Ti,pi (etat intern)
CC    OUT: ro0,u0,v0,w0,T0,pi (etat extern)
CC    OUT: tcx,tcy,tcz, tnx,tny,tnz
#include  "FastS/BC/linearization_state.for"
CC    In : tnx,tny,tnz, c_ale, ventijk
CC    OUT: qen (vitesse entrainememnt)
#include  "FastS/BC/qn_ale.for"

C... End of Newton and building of the variables
          vng  = (wng(lit)+qen)
C...      Absolute velocity module
          vg   = vng*usd0n(lit)
          
          b    = 1. - (vng**2)*usd0n2(lit)/(2.*ha(li))
          b = max(0.2,b) !cutoff robustesse

          pg   = pa(li)*b**gam2
          rog  = gam2*pg/(ha(li)*b)
          
          Tg   = pg/(rog*rgp)

          !if(j*k.eq.1)write(*,'(a,4f16.6,2i4)')'verif',vg,pg,rog,Tg,j,k
          rop(l,1) = rog
          rop(l,2) = vg*d0x(li)
          rop(l,3) = vg*d0y(li)
          rop(l,4) = vg*d0z(li)
          rop(l,5) = Tg
