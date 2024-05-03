            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent


            !dpdn=0  et T imposee (density deduit de T et p)
            rho_w = rop(ldjr,1)* rop(ldjr,5)/tw(li)
            rop(l,1) = 2*rho_w-rop(ldjr,1)

            u  =  rop(ldjr,2)
            v  =  rop(ldjr,3)
            w  =  rop(ldjr,4)

c            rop(l,1) = rop(ldjr,1)
            rop(l,2) = ventx*c_ale - u
            rop(l,3) = venty*c_ale - v
            rop(l,4) = ventz*c_ale - w

            !pour imposer la bonne temperature paroi pour les terme
            !visqueux, il faut que 0.5*(Tfictif+Treel)=tw
            ! donc:
            rop(l,5) = 2*tw(li)-rop(ldjr,5)

            !rop(l,5) = tw(li)
