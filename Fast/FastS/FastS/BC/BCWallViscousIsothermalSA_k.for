            ldjr = inddm( i, j, kref - k)
            li   = indbci(i, j)
#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for" 
             rop(l,6) =-rop(ldjr,6)
