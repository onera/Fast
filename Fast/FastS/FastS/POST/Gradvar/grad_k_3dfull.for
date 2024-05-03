        l5 = l+inck
        l6 = l-inck

        tix = tk(lt,1)
        tiy = tk(lt,2)
        tiz = tk(lt,3)

        !!GradU
        u =(rop(l )+rop(l6 ))*c1- (rop(l5)+rop(l-inc2k))*c2

        dvardc( l  + V1 ) = dvardc( l  + V1 ) - u*tix
        dvardc( l0 + V1 ) = dvardc( l0 + V1 ) + u*tix
        dvardc( l  + V2 ) = dvardc( l  + V2 ) - u*tiy
        dvardc( l0 + V2 ) = dvardc( l0 + V2 ) + u*tiy
        dvardc( l  + V3 ) = dvardc( l  + V3 ) - u*tiz
        dvardc( l0 + V3 ) = dvardc( l0 + V3 ) + u*tiz

