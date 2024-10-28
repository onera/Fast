      ! on amorti les fluctu en dehors de la couche limite: si ijk > jamor defini dans muse
      amor  = 1
      iamor = ci_amor*i +cj_amor*j + ck_amor*k

      if (iamor.gt.jamor) then
        amor=abs(ijk_amor-iamor)/float(ijk_amor+1-jamor)
        amor=amor*amor
      endif

      roe     =rho_in(li)+ clund*(rop(lr,1)-rho_planlund(li))*amor
      roe_inv = 1./roe       

      rue = roe*( u_in(li) + clund*(rop(lr,2)- u_planlund(li))*amor )
      rve = roe*( v_in(li) + clund*(rop(lr,3)- v_planlund(li))*amor )
      rwe = roe*( w_in(li) + clund*(rop(lr,4)- w_planlund(li))*amor )

      ete = roe*t_in(li)*cv +0.5*( rue*rue+rve*rve+rwe*rwe)*roe_inv

#include  "FastS/BC/non_reflection.for"

      roinv = 1./svar1

      u  = u + r*svar2*roinv
      v  = v + r*svar3*roinv
      w  = w + r*svar4*roinv

      rop(l,1) =  svar1
      rop(l,2) =  u 
      rop(l,3) =  v 
      rop(l,4) =  w 
      rop(l,5) =  (qvar5*roinv - .5*(u*u+v*v+w*w))*cvinv

