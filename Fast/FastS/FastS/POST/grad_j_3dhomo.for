        l5 = l+incj
        l6 = l-incj

        tix = tj(lt,1)
        tiy = tj(lt,2)

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2j+V2))*c2

        dvardc( l  + Vdudx ) = dvardc( l  + Vdudx ) - u*tix
        dvardc( l0 + Vdudx ) = dvardc( l0 + Vdudx ) + u*tix
        dvardc( l  + Vdudy ) = dvardc( l  + Vdudy ) - u*tiy
        dvardc( l0 + Vdudy ) = dvardc( l0 + Vdudy ) + u*tiy

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2j+V3))*c2

        dvardc( l  + Vdvdx ) = dvardc( l  + Vdvdx ) - u*tix
        dvardc( l0 + Vdvdx ) = dvardc( l0 + Vdvdx ) + u*tix
        dvardc( l  + Vdvdy ) = dvardc( l  + Vdvdy ) - u*tiy
        dvardc( l0 + Vdvdy ) = dvardc( l0 + Vdvdy ) + u*tiy

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2j+V4))*c2

        dvardc( l  + Vdwdx ) = dvardc( l  + Vdwdx ) - u*tix
        dvardc( l0 + Vdwdx ) = dvardc( l0 + Vdwdx ) + u*tix
        dvardc( l  + Vdwdy ) = dvardc( l  + Vdwdy ) - u*tiy
        dvardc( l0 + Vdwdy ) = dvardc( l0 + Vdwdy ) + u*tiy


