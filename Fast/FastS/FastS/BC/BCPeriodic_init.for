      ! mise a zero si tableau drodm du RHS implicit
      if(lrhs.eq.1) then
        c1= 0.
      else !extrapolation ordre 0 si variable conservative ou primitive
        c1= 1.
      endif

