      Q(l) = dvardc(lg+ Vdudx)*dvardc(lg+ Vdudx)
      Vgrad= Vdvdx+dim_grad                            !!dvdy
      Q(l) = Q(l)+ dvardc(lg+ Vgrad)*dvardc(lg+ Vgrad)

       Q(l)=-0.5*Q(l)- dvardc(lg+ Vdudx+dim_grad)*dvardc(lg+ Vdvdx)   !dudy*dvdx

       !valeur interface green 2 fois trop grande(gradient auu carree
       !4 fois trop garnd
        Q(l)= Q(l)*volinv
