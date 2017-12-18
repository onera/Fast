#include  "FastS/Compute/mulam.for" 
#include  "FastS/Compute/SA/chi.for"
          fvv1 = fv1(chi)
#include  "FastS/Compute/SA/xmut.for" 
          ad1      =dlng(l)
          adcut    =max(ad1 , 1.e-27)
          ra       =(xmut(l)/rop(l ,1))
     1                  /(SA_CKARM*SA_CKARM*adcut*adcut*auijuij**0.5)
          variable2=512.*ra*ra*ra
          variable2=min(variable2,15.)
          fa       =1. - tanh(variable2)
          testfa   =0.5+ sign(.5, 0.8-fa)
          adelta1  =testfa*sph2 +(1.-testfa)*adelta1                    ! echelle caractéristique maillage

          !choix distance
          adtild1  =ad1-fa*max(ad1-SA_RCTEDES*adelta1,0.)

          dist     =max(adtild1,1.e-27)

          !CALCUL DU TERME DE PRODUCTION
          f2       = fv2(chi,fvv1)

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

          tsource   = (prod+anvisc-destruc)
