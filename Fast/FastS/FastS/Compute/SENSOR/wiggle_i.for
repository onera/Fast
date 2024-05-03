#include  "FastS/Compute/SENSOR/wiggle.for"

         test=min(f1,f2,f3,f4,f5)
         if(test.lt.souszero) then
            wig(l +v1)   = 1.
         else
            wig(l +v1)   = 0.
         endif
         !wig(l +v1) = param_real(37)
