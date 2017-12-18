        l5 = l+inck
        l6 = l-inck

        tiz = tk(lt,1)

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2k+V2))*c2

        dvardc( ls + Vdudz ) = dvardc( ls + Vdudz ) + u*tiz

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2k+V3))*c2

        dvardc( ls + Vdvdz ) = dvardc( ls + Vdvdz ) + u*tiz

        !!GradW
        u =(rop(l +V4)+rop(l6 +V4))*c1- (rop(l5+V4)+rop(l-inc2k+V4))*c2

        dvardc( ls + Vdwdz ) = dvardc( ls + Vdwdz ) + u*tiz


