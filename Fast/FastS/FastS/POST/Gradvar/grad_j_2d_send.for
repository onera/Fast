        l5 = l+incj
        l6 = l-incj

        tix = tj(lt,1)
        tiy = tj(lt,2)

        !!GradU
        u =(rop(l)+rop(l6))*c1- (rop(l5)+rop(l-inc2j))*c2

        dvardc( ls + V1 ) = dvardc( ls + V1 ) + u*tix
        dvardc( ls + V2 ) = dvardc( ls + V2 ) + u*tiy

