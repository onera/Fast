            ldjr = inddm( i, jref - j, k)
            li   = indbci(i,  k )
#include     "FastS/BC/BCWallViscousIsothermal_firstrank.for"
            rop(l,6) =-rop(ldjr,6)
