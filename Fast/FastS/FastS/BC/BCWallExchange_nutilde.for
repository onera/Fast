      yplus = rho_w*utau*xmut(ltg , 2)/mu_w
      !van driest pour nut
      expy     = 1.-exp(-yplus/19.)! ranges 17 a 26
      !nutcible = abs( kappa*aa*utau*expy**2)!negatif si pt ibc interieur aux corps
      nutcible = abs( kappa*rop(ltg , 1)*xmut(ltg , 2)*utau*expy**2)!negatif si pt ibc interieur aux corps
      !equation 4eme degre
      nu = (coesut*sqrt(rop(ltg , 5))*rop(ltg , 5)/(rop(ltg , 5)+cmus1))
     &    /rop(ltg , 1)
      a = nutcible/nu
      b = a*cv1*cv1*cv1

      !changement de variable 4eme degre
      p = -3*a*a/8.
      q = -a*a*a/8.
      r = -3*a*a*a*a/256.-b
      !a_adim = a/(mu_vec[noind]/ro_vec[noind])

      !equation 3eme degre
      ap = 8.
      bp = -4*p
      cp = -8*r
      dp = 4*p*r-q*q

      !racines du 3eme degre
      delta = 18.*ap*bp*cp*dp-4.*bp*bp*bp*dp+bp*bp*cp*cp
     &        -4.*ap*cp*cp*cp-27.*ap*ap*dp*dp
      delta0 = bp*bp-3.*ap*cp
      delta1 = 2.*bp*bp*bp-9.*ap*bp*cp+27.*ap*ap*dp

      superdelta = -27.*ap*ap*delta

      y1 = -1.
      y2 = -1.
      if (abs(delta).lt.tol.and.abs(delta0).lt.tol) then 
        y1 = -bp/(3*ap)
      elseif (abs(delta).lt.tol) then
        y1 = (9*ap*dp - bp*cp)/(2*delta0)
        y2 = (4*ap*bp*cp-9*ap*ap*dp-bp*bp*bp)/(ap*delta0)
      else
         !version super delta
         c1 = -1.
         c2 = -1.
         if (superdelta.ge.0.) then
           root = sqrt(superdelta)
           if (delta1-root.ge.0.) then
             c1 = ((delta1 -root)*0.5)**(1./3.)
           else
             c1 = -((root -delta1)*0.5)**(1./3.)
           endif
           if (delta1+root.ge.0.) then
             c2 = ((delta1 +root)*0.5)**(1./3.)
           else
             c2 = -((-root-delta1)*0.5)**(1./3.)
           endif
         endif
         y1 = -1./(3*ap)*(bp + c1 + delta0/c1 )
         y2 = -1./(3*ap)*(bp + c2 + delta0/c2 )
         !printf("racine y1=%g y2=%g\n", y1, y2);
      endif
 

      !racine de l equation du 4eme degre
      c1 = 2*y1-p
      c2 = 2*y2-p
      !printf("c1 > 0 = %g, c2 > 0 = %g\n",c1,c2);

      z1 = -123456.
      if (c1.ge.tol) then
        p1 = -2*y1-p+2*q/(sqrt(c1))
        p2 = -2*y1-p-2*q/(sqrt(c1))
        !printf("1. p1=%g p2=%g\n", p1,p2);
        if (p1.ge.tol) then
          z1 = 0.5*(sqrt(c1)+sqrt(p1))
        elseif (p2.ge.tol) then 
          z1 = 0.5*(sqrt(c1)+sqrt(p2))
        endif
      endif

      if (c2.ge.tol.and.z1.eq.-123456) then
        p1 = -2*y2-p+2*q/(sqrt(c2))
        p2 = -2*y2-p-2*q/(sqrt(c2))
        if (p1.ge.tol) then
          z1 = 0.5*(sqrt(c2)+sqrt(abs(p1)))
        elseif (p2.ge.tol) then
          z1 = 0.5*(sqrt(c2)+sqrt(abs(p2)))
        endif
      endif
      if (c1.le.tol.and.z1.eq.-123456) then
        b0 = y1*y1-r
        if (b0.ge.tol) then
          p1 = -2*y1-p+4.*sqrt(b0)
          p2 = -2*y1-p-4.*sqrt(b0)
          if (p1.ge.tol) then 
             z1 = 0.5*(sqrt(abs(c1))+sqrt(abs(p1)))
          elseif (p2.ge.tol) then
             z1 = 0.5*(sqrt(abs(c1))+sqrt(abs(p2)))
          endif
        endif
      endif
      if (c2.le.tol.and.z1.eq.-123456) then
        b0 = y2*y2-r
        if (b0.ge.tol) then 
          p1 = -2*y2-p+4.*sqrt(b0)
          p2 = -2*y2-p-4.*sqrt(b0)
          if (p1.ge.tol) then
            z1 = 0.5*(sqrt(abs(c2))+sqrt(abs(p1)))
          elseif (p2 >= tol) then
            z1 = 0.5*(sqrt(abs(c2))+sqrt(abs(p2)))
          endif
        endif
      endif

      !nutile final
      rop(ltg, 6) = (z1 + a*0.25)*nu
