        l5 = l+incj
        l6 = l-incj

        !!GradU
        u =(rop(l )+rop(l6 ))*c1- (rop(l5)+rop(l-inc2j))*c2

        dvardc( l  + V2 ) = dvardc( l  + V2 ) - u*tjy
        dvardc( l0 + V2 ) = dvardc( l0 + V2 ) + u*tjy

