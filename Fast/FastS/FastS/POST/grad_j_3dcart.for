        l5 = l+incj
        l6 = l-incj

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2j+V2))*c2

        dvardc( l  + Vdudy ) = dvardc( l  + Vdudy ) - u*tjy
        dvardc( l0 + Vdudy ) = dvardc( l0 + Vdudy ) + u*tjy

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2j+V3))*c2

        dvardc( l  + Vdvdy ) = dvardc( l  + Vdvdy ) - u*tjy
        dvardc( l0 + Vdvdy ) = dvardc( l0 + Vdvdy ) + u*tjy

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2j+V4))*c2

        dvardc( l  + Vdwdy ) = dvardc( l  + Vdwdy ) - u*tjy
        dvardc( l0 + Vdwdy ) = dvardc( l0 + Vdwdy ) + u*tjy


