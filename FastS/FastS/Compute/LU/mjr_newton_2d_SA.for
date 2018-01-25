        ro_old    = rop(l,1)
        u_old     = rop(l,2)
        v_old     = rop(l,3)
        t_old     = rop(l,5)
        nu_old    = rop(l,6)

        rop_1(l,1) = rop(l,1) + drodm(l,1)
        r_1        = 1./rop_1(l,1)

        rop_1(l,2) = (ro_old*u_old + drodm(l,2))*r_1
        rop_1(l,3) = (ro_old*v_old + drodm(l,3))*r_1
        rop_1(l,4) = 0.

        roe_old        = ro_old*( param_real(CVINF)*t_old 
     &                           + 0.5*( u_old*u_old
     &                                  +v_old*v_old))

        rop_1(l,5)     = ( roe_old + drodm(l,5) )*r_1*cvinv
     &                  - cvinv2*( rop_1(l,2)*rop_1(l,2)
     &                            +rop_1(l,3)*rop_1(l,3))

        anulam = coesut*sqrt(rop_1(l,5)*temp01)/(1.+cmus1/rop_1(l,5))
        anulam = anulam*r_1*ratiom 

        t_old  = (ro_old*nu_old + drodm(l,6))*r_1

        t_old      = min(t_old,anulam)
        rop_1(l,6) = max(t_old,0.)
!        rop_1(l,6) = max(t_old,param_real(RONUTILDEINF))
