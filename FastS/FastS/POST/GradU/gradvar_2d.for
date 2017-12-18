        !u4 =(rop(l +Vvar)+rop(l4 +Vvar))
        !u3 =(rop(l +Vvar)+rop(l3 +Vvar))
        !u6 =(rop(l +Vvar)+rop(l6 +Vvar))
        !u5 =(rop(l +Vvar)+rop(l5 +Vvar))

        l5    = l+inci
        l6    = l-inci
        incdir=inc2i
#include   "FastS/POST/GradU/interp_muscl.for"
        dvardc(lg+ Vgrad           )= -u1*ti(lt+v1mtr) 
     &                                +u3*ti(lt+inci_mtr+v1mtr)

        dvardc(lg+ Vgrad+dim_grad  )= -u1*ti(lt+v2mtr) 
     &                                +u3*ti(lt+inci_mtr+v2mtr)

        l5    = l+incj
        l6    = l-incj
        incdir=inc2j
#include   "FastS/POST/GradU/interp_muscl.for"

        dvardc(lg+ Vgrad           )=  dvardc(lg+ Vgrad           )
     &                                -u1*tj(lt+v1mtr) 
     &                                +u3*tj(lt+incj_mtr+v1mtr)

        dvardc(lg+ Vgrad+dim_grad  )= dvardc(lg+ Vgrad+dim_grad  )
     &                                -u1*tj(lt+v2mtr) 
     &                                +u3*tj(lt+incj_mtr+v2mtr)
