#include   "FastS/POST/GradU/interp_muscl.for"

      dvardc(lg+ Vgrad+dim_grad  )=dvardc(lg+ Vgrad+dim_grad  )
     &                            - (u1 -u3)*tcy

