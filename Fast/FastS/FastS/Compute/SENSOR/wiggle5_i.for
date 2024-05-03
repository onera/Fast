#include  "FastS/Compute/SENSOR/wiggle5.for"

         test=min(f1,f2,f3,f4,f5)
         !test=min(f1,f2,f3,f5)
         if(test.lt.souszero) then
            wig(l +v1)   = 1.
         else if(test.lt.(1.0+souszero)) then 
            wig(l +v1)   = 0.6
         else if(test.lt.(3.0+souszero)) then 
            wig(l +v1)   = 0.4
         else 
            wig(l +v1)   = 0.05
         endif
