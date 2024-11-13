            ldjr = inddm( i, j, kr        )
            ldl  = inddm( i, j, kr +sens_int  )
            ldnp = inddm( i, j, kr +sens_int*2)
            ldnm = ldjr
 
            m    = ldjr  - exchange
            loo  = ldjr  - exchange*shift_loo
#include     "FastS/BC/BCWallModel.for"

