#include  "FastS/Compute/SENSOR_INIT/wiggle.for"

         test=min(f1,f2,f3,f4,f5)
         if(test.lt.souszero) then
            wig(l +v1)   = 1.
         else
            wig(l +v1)   = 0.
         endif
