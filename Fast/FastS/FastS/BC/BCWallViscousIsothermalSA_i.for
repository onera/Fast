            ldjr = inddm(  iref - i, j, k)
            li   = indbci(j, k)  
#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"
            rop(l,6) =-rop(ldjr,6)
            
