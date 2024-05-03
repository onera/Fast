        l5 = l+incj
        l6 = l-incj

        tix = tj(lt,1)
        tiy = tj(lt,2)

        !!GradU
        u =(rop(l)+rop(l6))*c1- (rop(l5)+rop(l-inc2j))*c2

        dvardc( l  + V1 ) = dvardc( l  + V1 ) - u*tix
        dvardc( l0 + V1 ) = dvardc( l0 + V1 ) + u*tix
        dvardc( l  + V2 ) = dvardc( l  + V2 ) - u*tiy
        dvardc( l0 + V2 ) = dvardc( l0 + V2 ) + u*tiy

