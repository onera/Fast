#if __DEBUG__
            ! encadre le depassement intentionnel
            if (c_ale > 0) then
#endif
            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent
#if __DEBUG__
            else
            ventx = 0.
            venty = 0.
            ventz = 0.
            endif
#endif
            qen = (ventx*tnx + venty*tny + ventz*tnz)*c_ale

            

