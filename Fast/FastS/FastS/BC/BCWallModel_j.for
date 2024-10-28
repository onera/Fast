            ldjr = inddm( i, jr        , k )
            ldl  = inddm( i, jr +sens  , k )
            ldnp = inddm( i, jr +sens*2, k )
            ldnm = ldjr
 
            m    = ldjr  - exchange
            loo  = ldjr  - exchange*shift_loo
#include     "FastS/BC/BCWallModel.for"

