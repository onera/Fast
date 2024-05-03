        l5 = l+inci
        l6 = l-inci

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2i+V2))*c2

        dvardc( l  + Vdudx ) = dvardc( l  + Vdudx ) - u*tix
        dvardc( l0 + Vdudx ) = dvardc( l0 + Vdudx ) + u*tix

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2i+V3))*c2

        dvardc( l  + Vdvdx ) = dvardc( l  + Vdvdx ) - u*tix
        dvardc( l0 + Vdvdx ) = dvardc( l0 + Vdvdx ) + u*tix

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2i+V4))*c2

        dvardc( l  + Vdwdx ) = dvardc( l  + Vdwdx ) - u*tix
        dvardc( l0 + Vdwdx ) = dvardc( l0 + Vdwdx ) + u*tix


