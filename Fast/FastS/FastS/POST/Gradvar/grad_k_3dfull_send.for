        l5 = l+inck
        l6 = l-inck

        tix = tk(lt,1)
        tiy = tk(lt,2)
        tiz = tk(lt,3)

        !!GradU
        u =(rop(l )+rop(l6 ))*c1- (rop(l5)+rop(l-inc2k))*c2

        dvardc( ls + V1 ) = dvardc( ls + V1 ) + u*tix
        dvardc( ls + V2 ) = dvardc( ls + V2 ) + u*tiy
        dvardc( ls + V3 ) = dvardc( ls + V3 ) + u*tiz
