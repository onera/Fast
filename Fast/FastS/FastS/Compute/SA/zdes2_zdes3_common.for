          stild    = rot + (anutild*f2)/(SA_CKARM*SA_CKARM*dist*dist) ! \hat{2}=\Omega + \hat{\nu}*fv_2/(\kappa^2*d^2)
          
          !!Correction SA-R
          stild    = stild + SWITCH_SA_ROT_CORR
     &               *SA_CROT*min(0.00000000000000000001,St-rot)      ! \hat{S}+C_{rot}*min(0,\{hat}-\Onera)

          prod     = rop(l,1)*SA_CB1*stild*anutild                    ! 
       
          !CALCUL DU TERME DE DESTRUCTION
          stild     = max(stild,0.00000000000000000001)             !\hat{S}
          r         = anutild/(stild*SA_CKARM*SA_CKARM*dist*dist)   !\hat{v]/(\hat{S}*\kapa^2*d^2)
          r         = min(r,10.)                                    !min(\hat{v]/(\hat{S}*\kapa^2*d^2),10)
          r5        = (r*r*r*r*r-1.)                                !r^5-1
          g         = r *(1.+SA_CW2*r5)                             !g=r+c_{w2}*(r^6-1)
          
          !!Correction SA-LRE
          SA_CW2_LRE = SA_CW4 + SA_CW5 / ( (1.+ chi/40.)**2)        !ok
          g  = SWITCH_SA_LOW_RE*r*(1.+SA_CW2_LRE*r5) + 
     &               (1.-SWITCH_SA_LOW_RE)*g                        !ok
