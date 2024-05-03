      ! Wall function for the temperature based on
      ! "An adjustment to the standard temperature wall function for
      ! CFD modeling of indoor convective heat transfer"
      ! Building and environment. Zhang et al (2013)

      uplus = unorm/utau

      dist = xmut(1+param_int(NDIMDX)) !wall distance
      yplus = rop(m+v1)*dist*utau/xmut(m)

      pr = cp*xmut(m)/cond
      prtur = 0.9
      ratio = pr/prtur

      ! Pee function (Jayatilleke 1969)
      pf    = 9.24*(ratio**0.75-1.0)*(1.0+0.28*exp(-0.007*ratio))

      ! Blending function dt
      dt    = - (0.01*(pr*yplus)**4)/(1+5*yplus*pr**3)

      Tplus = (pr*yplus)*exp(dt)+prtur*(uplus+pf)*exp(1/dt)

      qwall = rop(m+v1)*cp*utau*(Twall-rop(m+v5))/Tplus
