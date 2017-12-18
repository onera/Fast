            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent
            
            !write(*,*) ventx, venty, ventz, j,k
          
            qen = (ventx*tnx + venty*tny + ventz*tnz)*c_ale

            

