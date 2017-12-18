        l5 = l+inci
        l6 = l-inci

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2i+V2))*c2

        dvardc( ls + Vdudx ) = dvardc( ls + Vdudx ) + u*tix

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2i+V3))*c2

        dvardc( ls + Vdvdx ) = dvardc( ls + Vdvdx ) + u*tix

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2i+V4))*c2

        dvardc( ls + Vdwdx ) = dvardc( ls + Vdwdx ) + u*tix


