                  ldjr   = inddm( i, jr , k )
                  lsample= ldjr + (sample-1)*incj
                  l      = inddm( i, j+2*incj , k )
#include          "FastS/BC/BCWallExchange.for"

                  !remplissage ghost pour visu
                  !lghost   = inddm( i , j-2  , k )
                  lghost   = inddm( i , j  , k )
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)

