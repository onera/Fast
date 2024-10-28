      sampling = param_int(WM_SAMPLING)

      exchange = 1 - sampling  ! index correction for identifying the correct cell

      sens     = 1.
      sens_int = 1.
      shift = 0
      if(mod(idir,2).eq.0) then
        sens      = -1.
        sens_int  = -1
        shift = 1 
        exchange =  -exchange
      endif

      shift =incijk*shift
      shiftvent = param_int(NIJK_VENT + idir/2)* exchange

      if(idir.le.2) then
        incijk=1
      elseif(idir.le.4) then
        incijk=param_int(NIJK)
      else
        incijk=param_int(NIJK+1)*param_int(NIJK)
      endif

      exchange =incijk*exchange
      shift_loo= 2
