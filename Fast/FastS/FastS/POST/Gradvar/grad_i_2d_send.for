        l5 = l+inci
        l6 = l-inci

        tix = ti(lt,1)
        tiy = ti(lt,2)

        !!GradU
        u =(rop(l)+rop(l6))*c1- (rop(l5)+rop(l-inc2i))*c2

        dvardc( ls + V1 ) = dvardc( ls + V1 ) + u*tix
        dvardc( ls + V2 ) = dvardc( ls + V2 ) + u*tiy

