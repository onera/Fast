#include  "FastS/Compute/mulam.for" 
#include  "FastS/Compute/SA/chi_diff.for" 

          !CALCUL DU TERME DE PRODUCTION
          f1       = fv1(chi)
          f2       = fv2(chi,f1)
          
          dist     = max(dlng(l), 1.e-27)

          !Version F Renac
          rot = rot+ sign(1.e-14, rot)
          stild    = (anutild*f2)/(SA_CKARM*SA_CKARM*dist*dist)
          !Oliver modification
          if(stild+SA_CV2*rot.ge.0.)then
             stild = rot + stild
          else
             stild = rot*(1. + (SA_CV2*SA_CV2*rot + SA_CV3*stild)
     &                        /(SA_CV3-2.*SA_CV2*rot-stild)
     &                   )
          endif

          stild    = rot + (anutild*f2)/(SA_CKARM*SA_CKARM*dist*dist)

          !!Correction SA-R (desactivee car St pas calcule en diff)
          !stild = stild + param_int(SA_ROT_CORR)
     &    !        *SA_CROT*min(0.00000000000000000001,St-rot) 

          prod     = rop(l,1)*SA_CB1*stild*anutild
       
          !CALCUL DU TERME DE DESTRUCTION
          !Version F Renac: pas de limitation stild
          r         = anutild/(stild*SA_CKARM*SA_CKARM*dist*dist)
          r         = min(r,10.)
          r5        = (r*r*r*r*r-1.)
          g         = r *(1.+SA_CW2*r5)
          !!Correction SA-LRE
          SA_CW2_LRE = SA_CW4 + SA_CW5 / ( (1.+ chi/40.)**2) 
          g  = param_int(SA_LOW_RE)*r*(1.+SA_CW2_LRE*r5) +
     &               (1-param_int(SA_LOW_RE))*g
          fwg       = fw(g)

          destruc   =rop(l,1)*cw1*fwg*(anutild/dist)*(anutild/dist)

          !CALCUL DU TERME SOURCE GLOBAL

          tsource   = (prod+anvisc-destruc)
