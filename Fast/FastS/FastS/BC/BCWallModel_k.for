            ldjr = inddm( i, j, kr        )
            ldl  = inddm( i, j, kr +sens  )
            ldnp = inddm( i, j, kr +sens*2)
            ldnm = ldjr
 
            m    = ldjr  - exchange
            loo  = ldjr  - exchange*shift_loo
#include     "FastS/BC/BCWallModel.for"

