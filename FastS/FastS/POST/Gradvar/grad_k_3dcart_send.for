        l5 = l+inck
        l6 = l-inck

        !!GradU
        u =(rop(l )+rop(l6 ))*c1- (rop(l5)+rop(l-inc2k))*c2

        dvardc( ls + V3 ) = dvardc( ls + V3 ) + u*tkz

