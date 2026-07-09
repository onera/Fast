#include "FastS/Compute/SENSOR/wiggle.for"

         test=min(f1,f2,f3,f4,f5)
         !filtrage temporel: wig(n+1) = c2 wig(n) + c1*wig(n+1)
         if(test.lt.souszero) then
            !wig(l +v3)   = 1.
            wig(l +v3)   = c2*wig(l +v3) + c1
         else
            !wig(l +v3)   = 0.
            wig(l +v3)   = c2*wig(l +v3)
         endif
