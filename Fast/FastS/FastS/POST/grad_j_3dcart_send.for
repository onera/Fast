        l5 = l+incj
        l6 = l-incj

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2j+V2))*c2

        dvardc( ls + Vdudy ) = dvardc( ls + Vdudy ) + u*tjy

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2j+V3))*c2

        dvardc( ls + Vdvdy ) = dvardc( ls + Vdvdy ) + u*tjy

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2j+V4))*c2

        dvardc( ls + Vdwdy ) = dvardc( ls + Vdwdy ) + u*tjy


