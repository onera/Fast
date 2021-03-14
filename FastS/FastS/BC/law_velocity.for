      ! ilaw = 1:
      ! Wall function for the velocity proposed in
      ! "Explicit expression for the smooth wall velocity
      ! distribution in a turbulent boundary layer"
      ! AIAA Journal. Musker (1979)

      ! ilaw =/ 1:
      ! Power law of the wall (Werner and Wengle)

      dist = xmut(m+param_int(NDIMDX)) !wall distance

      utau = ( ((unorm**7)*(xmut(m)/rop(m+v1)))/((8.3**7)*dist) )**0.125

      ilaw = 2
      if (ilaw.eq.1) then

         ! Newton method to determine utau using Musker's wall law
         iter = 0 
         f = 1.0
         do while (abs(f).gt.1e-6.and.iter.le.100)
            iter = iter + 1

            aa = rop(m+v1)*dist/xmut(m)
            yplus = aa*utau
            l1 = (yplus+10.6)**9.6
            l2 = yplus*(yplus-8.15) + 86
            l3 = (2*yplus-8.15)/16.7

            f = 5.424*atan(l3) + log10(l1/(l2*l2)) -
     &          unorm/utau - 3.52

            tp =  aa*(9.6*((yplus + 10.6)**8.6)  
     &          - l1/l2*( 4*yplus - 16.30 ) )

            fp = 0.649580838323*aa/(1 + l3*l3 )  +  tp/(l1*2.302585093)
     &          + unorm/(utau*utau)

            utau = utau - f/fp
         end do
      end if
      tauw = r*utau**2
