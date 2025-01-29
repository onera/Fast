                flux(l0             ) = flu1*sens
                flux(l0 +size_flux  ) = flu2*sens
                flux(l0 +size_flux*2) = flu3*sens
                flux(l0 +size_flux*3) = flu4*sens
                flux(l0 +size_flux*4) = flu5*sens
                !write(*,*)'adr', l0, l0 +size_flux*4
                !flux(l0             ) = sens
                !flux(l0 +size_flux  ) = 2*sens
                !flux(l0 +size_flux*2) = 3*sens
                !flux(l0 +size_flux*3) = 4*sens
                !flux(l0 +size_flux*4) = 5*sens
