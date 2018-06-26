        ro_old     = rop(l,1)
        u_old      = rop(l,2)
        v_old      = rop(l,3)
        w_old      = rop(l,4)
        t_old      = rop(l,5)

        rop_1(l,1) = rop(l,1) + drodm_out(l,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm_out(l,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm_out(l,3))*r_1
        rop_1(l,4) = (ro_old*w_old + drodm_out(l,4))*r_1

        roe_old        = ro_old*( param_real(CVINF)*t_old 
     &                           + 0.5*( u_old*u_old
     &                                  +v_old*v_old
     &                                  +w_old*w_old))

        rop_1(l,5)     = ( roe_old +  drodm_out(l,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3)
     &                            +rop_1(l,4)*rop_1(l,4))
