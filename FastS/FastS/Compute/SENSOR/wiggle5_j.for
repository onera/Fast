#include  "FastS/Compute/SENSOR/wiggle5.for"

         test=min(f1,f2,f3,f4,f5)
         if(test.lt.souszero) then
            wig(l +v2)   = 1.
         else if(test.lt.(1.0+souszero)) then 
            wig(l +v2)   = 0.6
         else if(test.lt.(2.0+souszero)) then 
            wig(l +v2)   = 0.5
         else if(test.lt.(3.0+souszero)) then 
            wig(l +v2)   = 0.4
         else if(test.lt.(4.0+souszero)) then 
            wig(l +v2)   = 0.3
         else 
            wig(l +v2)   = 0.05
         endif

