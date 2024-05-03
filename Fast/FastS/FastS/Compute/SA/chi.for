        anulam  = rop(l ,1)/amulam

        anutild = max(rop(l,6)      , 1.e-27)
        chi     = max(anutild*anulam, 1.e-7 ) 

        !Version F Renac
        !anutild = rop(l,6)
        !chi = anutild*anulam
        !if(abs(chi).lt.1.5)THEN
        !  chi = 0.05*LOG(1.+exp(20.*chi))
        !elseif(chi.lt. 0.)THEN
        !  chi = 0.
        !endif
