                  ldjr   = inddm( i ,  j , kr )
                  lsample= ldjr + (sample-1)*inck
                  l      = inddm( i, j, k+2*inck )
#include          "FastS/BC/BCWallExchangeSA.for"

                  !remplissage ghost pour visu
                  !lghost   = inddm( i , j-2  , k )
                  lghost   = inddm( i , j  , k )
                  rop(lghost,1) = rop(ldjr,1) 
                  rop(lghost,2) = rop(ldjr,2)
                  rop(lghost,3) = rop(ldjr,3)
                  rop(lghost,4) = rop(ldjr,4) 
                  rop(lghost,5) = rop(ldjr,5)
                  rop(lghost,6) = rop(ldjr,6)

