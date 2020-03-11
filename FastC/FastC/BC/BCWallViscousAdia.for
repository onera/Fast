            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent

           !on calcule l'inverse de la surface de la facette

            u  =  rop(il,2)
            v  =  rop(il,3)
            w  =  rop(il,4)

            rop(ir,1) = rop(il,1)
            rop(ir,2) = ventx*c_ale - u
            rop(ir,3) = venty*c_ale - v
            rop(ir,4) = ventz*c_ale - w

            rop(ir,5) = rop(il,5)
