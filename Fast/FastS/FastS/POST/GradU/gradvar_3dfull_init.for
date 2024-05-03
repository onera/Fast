#include   "FastS/POST/GradU/interp_muscl.for"

        dvardc(lg+ Vgrad           )= - u1*tcx +u3*tcx1
        dvardc(lg+ Vgrad+dim_grad  )= - u1*tcy +u3*tcy1
        dvardc(lg+ Vgrad+dim_grad*2)= - u1*tcz +u3*tcz1
