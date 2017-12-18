        l5 = l+inci
        l6 = l-inci

        tix = ti(lt,1)
        tiy = ti(lt,2)

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2i+V2))*c2

        dvardc( ls + Vdudx ) = dvardc( ls + Vdudx ) + u*tix
        dvardc( ls + Vdudy ) = dvardc( ls + Vdudy ) + u*tiy

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2i+V3))*c2

        dvardc( ls + Vdvdx ) = dvardc( ls + Vdvdx ) + u*tix
        dvardc( ls + Vdvdy ) = dvardc( ls + Vdvdy ) + u*tiy


