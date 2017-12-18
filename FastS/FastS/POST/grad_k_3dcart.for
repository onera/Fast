        l5 = l+inck
        l6 = l-inck

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2k+V2))*c2

        dvardc( l  + Vdudz ) = dvardc( l  + Vdudz ) - u*tkz
        dvardc( l0 + Vdudz ) = dvardc( l0 + Vdudz ) + u*tkz

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2k+V3))*c2

        dvardc( l  + Vdvdz ) = dvardc( l  + Vdvdz ) - u*tkz
        dvardc( l0 + Vdvdz ) = dvardc( l0 + Vdvdz ) + u*tkz

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2k+V4))*c2

        dvardc( l  + Vdwdz ) = dvardc( l  + Vdwdz ) - u*tkz
        dvardc( l0 + Vdwdz ) = dvardc( l0 + Vdwdz ) + u*tkz


