#include  "FastS/Compute/mulam.for" 
#include  "FastS/Compute/SA/chi.for"
          !choix distance
          ad1     = max( dlng(l), 1.e-27)

          testfa  =0.5+ sign(.5,SA_RCTEDES*adelta1 - ad1 )
          fvv1    = (1.-testfa)+testfa*fv1(chi)
#include  "FastS/Compute/SA/xmut.for"

          adtild1 = min( ad1    , SA_RCTEDES*adelta1)
          dist    = max( adtild1, 1.e-27)

          !CALCUL DU TERME DE PRODUCTION
          f2       = fv2(chi,fvv1)*testfa

          stild    = rot + (anutild*f2)/(SA_CKARM*SA_CKARM*dist*dist)

          prod     = rop(l,1)*SA_CB1*stild*anutild

       
          !CALCUL DU TERME DE DESTRUCTION
          stild     = max(stild,0.00000000000000000001)
          r         = anutild/(stild*SA_CKARM*SA_CKARM*dist*dist)
          r         = min(r,10.)
          g         = r *(1.+  SA_CW2*(r*r*r*r*r-1.) )
          ! verifier formule
          !fwg       = fw(g)*(1.-testfa) + testfa
          fwg       = (1.-testfa) + testfa*fw(g)

          destruc   =rop(l,1)*cw1*fwg*(anutild/dist)*(anutild/dist)

C       if(ndom.eq.1.and.j.eq.80.and.l-lij+ind_loop(1).eq.400.and.k.eq.1)
C     & write(*,'(5f20.12)') destruc, testfa, dist, adelta1, ad1
 
          !CALCUL DU TERME SOURCE GLOBAL

          tsource   = (prod+anvisc-destruc)
