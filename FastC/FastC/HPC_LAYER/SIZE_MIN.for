      !Calcul de la taille minimal 1D du bloc thread 
      if(param_int(IFLOW).eq.4) then !lbm
         lmin = 1
      elseif(param_int(ITYPCP).eq.2) then ! ns explicit
         lmin = 4
      else
         lmin =10
      endif
