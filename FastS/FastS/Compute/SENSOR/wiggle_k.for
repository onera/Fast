#include "FastS/Compute/SENSOR/wiggle.for"

         test=min(f1,f2,f3,f4,f5)
         if(test.lt.souszero) then
            wig(l +v3)   = 1.
         else
            wig(l +v3)   = 0.
         endif
