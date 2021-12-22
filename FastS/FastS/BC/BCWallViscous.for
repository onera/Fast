            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent

            u  =  rop(ldjr,2)
            v  =  rop(ldjr,3)
            w  =  rop(ldjr,4)

            rop(l,1) = rop(ldjr,1)
            rop(l,2) = 2.*ventx*c_ale - u
            rop(l,3) = 2.*venty*c_ale - v
            rop(l,4) = 2.*ventz*c_ale - w

            rop(l,5) = rop(ldjr,5)
