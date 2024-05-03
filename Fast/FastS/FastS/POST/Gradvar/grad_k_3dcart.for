        l5 = l+inck
        l6 = l-inck

        !!GradU
        u =(rop(l )+rop(l6 ))*c1- (rop(l5)+rop(l-inc2k))*c2

        dvardc( l  + V3 ) = dvardc( l  + V3 ) - u*tkz
        dvardc( l0 + V3 ) = dvardc( l0 + V3 ) + u*tkz

