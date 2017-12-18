      Q(l) = dvardc(lg+ Vdudx)*dvardc(lg+ Vdudx)
      Vgrad= Vdvdx+dim_grad                            !!dvdy
      Q(l) = Q(l)+ dvardc(lg+ Vgrad)*dvardc(lg+ Vgrad)

      Vgrad = Vdwdx+2*dim_grad                          !!dwdz
      Q(l) = Q(l)+ dvardc(lg+ Vgrad)*dvardc(lg+ Vgrad)
                
       Q(l)=-0.5*Q(l)
     &     - dvardc(lg+ Vdudx+dim_grad)*dvardc(lg+ Vdvdx)           !dudy*dvdx
     &     - dvardc(lg+ Vdwdx+dim_grad)*dvardc(lg+ Vdvdx+2*dim_grad)!dwdy*dvdz
     &     - dvardc(lg+ Vdwdx         )*dvardc(lg+ Vdudx+2*dim_grad)!dwdx*dudz

       !valeur interface green 2 fois trop grande(gradient auu carree
       !4 fois trop garnd
        Q(l)= Q(l)*volinv

