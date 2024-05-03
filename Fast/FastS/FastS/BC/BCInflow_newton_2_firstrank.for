CC    IN : rop
CC    OUT: roi,ui,vi,wi,Ti,pi (etat intern)
CC    OUT: ro0,u0,v0,w0,T0,pi (etat extern)
CC    OUT: tcx,tcy,tcz, tnx,tny,tnz
#include  "FastS/BC/linearization_state.for"
CC    In : tnx,tny,tnz, c_ale, ventijk
CC    OUT: qen (vitesse entrainememnt)
#include  "FastS/BC/qn_ale.for"
C... Newton for the relative normal velocity
          
          b = 1. - ((wng(lit)+qen)**2)*usd0n2(lit)/(2.*ha(li))

          b = max(0.2,b) !cutoff robustesse

          ! p = Pi (1 -U^2/(2CpTi))^(gamma/(gamma-1))
          pg  = pa(li)*b**gam2
 
c      nan = isnan(pg)
c      if(nan)   write(*,*)'fuck Nan inflow_newton',pa(li),b,gam2

          rog = gam2*pg/(ha(li)*b)

          f    = pg + roc0(lit)*wng(lit) - ri(lit)
          df   = roc0(lit) - rog*(wng(lit)+qen)*usd0n2(lit)
 
          dwng = -f/df 
          wng(lit) = wng(lit) + dwng
          
          residug = max(residug, ABS(dwng/wng(lit)) )

