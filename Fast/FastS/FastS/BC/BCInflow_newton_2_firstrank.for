C... Newton for the relative normal velocity
          
          b = 1. - ((wng+qen)**2)*usd0n2/(2.*ha(li))

          b = max(0.2,b) !cutoff robustesse

          ! p = Pi (1 -U^2/(2CpTi))^(gamma/(gamma-1))
          pg  = pa(li)*b**gam2
 
c      nan = isnan(pg)
c      if(nan)   write(*,*)'fuck Nan inflow_newton',pa(li),b,gam2

          rog = gam2*pg/(ha(li)*b)

          f    = pg + roc0*wng - ri
          df   = roc0 - rog*(wng+qen)*usd0n2
 
          dwng = -f/df 
          wng = wng + dwng
          
          residug = max(residug, ABS(dwng/wng) )

