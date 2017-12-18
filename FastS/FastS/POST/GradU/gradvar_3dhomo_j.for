#include   "FastS/POST/GradU/interp_muscl.for"

      dvardc(lg+ Vgrad           )=dvardc(lg+ Vgrad) 
     &                            - u1*tcx +u3*tcx1
      dvardc(lg+ Vgrad+dim_grad  )=dvardc(lg+ Vgrad+dim_grad  )
     &                            - u1*tcy +u3*tcy1

