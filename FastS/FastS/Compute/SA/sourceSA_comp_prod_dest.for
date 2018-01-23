#include  "FastS/Compute/mulam.for" 
#include  "FastS/Compute/SA/chi.for" 

          !CALCUL DU TERME DE PRODUCTION
          f1       = fv1(chi)
          f2       = fv2(chi,f1)
          
          dist     = max(dlng(l), 1.e-27)
          stild    = rot + (anutild*f2)/(SA_CKARM*SA_CKARM*dist*dist)

          prod     = rop(l,1)*SA_CB1*stild*anutild
       
          !CALCUL DU TERME DE DESTRUCTION
          stild     = max(stild,0.00000000000000000001)
          r         = anutild/(stild*SA_CKARM*SA_CKARM*dist*dist)
          r         = min(r,10.)
          g         = r *(1.+  SA_CW2*(r*r*r*r*r-1.) )
          fwg       = fw(g)

          destruc   =rop(l,1)*cw1*fwg*(anutild/dist)*(anutild/dist)

          !CALCUL DU TERME SOURCE GLOBAL

          tsource   = (prod+anvisc-destruc)+cc
