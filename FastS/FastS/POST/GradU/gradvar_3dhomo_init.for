#include   "FastS/POST/GradU/interp_muscl.for"

        dvardc(lg+ Vgrad           )= 0.
        dvardc(lg+ Vgrad+dim_grad  )= 0.
        dvardc(lg+ Vgrad+dim_grad*2)= - u1*tcz +u3*tcz1
