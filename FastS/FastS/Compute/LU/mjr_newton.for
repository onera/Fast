        ro_old     = rop(l,1)
        u_old      = rop(l,2)
        v_old      = rop(l,3)
        w_old      = rop(l,4)
        t_old      = rop(l,5)

        rop_1(l,1) = rop(l,1) + drodm_out(ls,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm_out(ls,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm_out(ls,3))*r_1
        rop_1(l,4) = (ro_old*w_old + drodm_out(ls,4))*r_1

        roe_old        = ro_old*( param_real(CVINF)*t_old 
     &                           + 0.5*( u_old*u_old
     &                                  +v_old*v_old
     &                                  +w_old*w_old))

        rop_1(l,5)     = ( roe_old +  drodm_out(ls,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3)
     &                            +rop_1(l,4)*rop_1(l,4))


c         rop_1(l,1) = max(0.25,rop_1(l,1))
c         rop_1(l,1) = min(3.9,rop_1(l,1))
c       
c         rop_1(l,2) = max(-500.0,rop_1(l,2) )
c         rop_1(l,2) = min( 600.0,rop_1(l,2) )
c
c         rop_1(l,3) = max(-650.0,rop_1(l,3) )
c         rop_1(l,3) = min( 650.0,rop_1(l,3) ) 

c         rop_1(l,4) = max(-600.0,rop_1(l,4) )
c         rop_1(l,4) = min( 520.0,rop_1(l,4) )

c         rop_1(l,5) = max( 100.0,rop_1(l,5) )
c         rop_1(l,5) = min( 600.0,rop_1(l,5) )
