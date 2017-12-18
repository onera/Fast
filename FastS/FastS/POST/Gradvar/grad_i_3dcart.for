        l5 = l+inci
        l6 = l-inci

        !!GradU
        u =(rop(l )+rop(l6 ))*c1- (rop(l5)+rop(l-inc2i))*c2

        dvardc( l  + V1 ) = dvardc( l  + V1 ) - u*tix
        dvardc( l0 + V1 ) = dvardc( l0 + V1 ) + u*tix


