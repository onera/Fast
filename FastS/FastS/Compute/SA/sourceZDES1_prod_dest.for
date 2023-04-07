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

          !!Correction SA-R
          St = 0.
          stild    = stild + SWITCH_SA_ROT_CORR
     &               *SA_CROT*min(0.00000000000000000001,St-rot) 

          prod     = rop(l,1)*SA_CB1*stild*anutild
       
          !CALCUL DU TERME DE DESTRUCTION
          stild     = max(stild,0.00000000000000000001)
          r         = anutild/(stild*SA_CKARM*SA_CKARM*dist*dist)
          r         = min(r,10.)
          r5        = (r*r*r*r*r-1.)
          g         = r *(1.+SA_CW2*r5)
          !!Correction SA-LRE
          SA_CW2_LRE = SA_CW4 + SA_CW5 / ( (1.+ chi/40.)**2) 
          g  = SWITCH_SA_LOW_RE*r*(1.+SA_CW2_LRE*r5) +
     &               (1.-SWITCH_SA_LOW_RE)*g
          ! verifier formule
          !fwg       = fw(g)*(1.-testfa) + testfa
          fwg       = (1.-testfa) + testfa*fw(g)

          destruc   =rop(l,1)*cw1*fwg*(anutild/dist)*(anutild/dist)

C       if(ndom.eq.1.and.j.eq.80.and.l-lij+ind_loop(1).eq.400.and.k.eq.1)
C     & write(*,'(5f20.12)') destruc, testfa, dist, adelta1, ad1
 
          !CALCUL DU TERME SOURCE GLOBAL

          tsource   = (prod+anvisc-destruc)
