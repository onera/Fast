      ! on amorti les fluctu en dehors de la couche limite: si ijk > jamor defini dans muse
      amor  = 1
      iamor = ci_amor*i +cj_amor*j + ck_amor*k

      if (iamor.gt.jamor) then
        amor=abs(ijk_amor-iamor)/float(ijk_amor+1-jamor)
        amor=amor*amor
      endif

      roe     =ro_inflow(li,1)+ clund*(rop(lr,1)-ro_planlund(li,1))*amor
      roe_inv = 1./roe       

      rue = roe*(   ro_inflow(li,2)
     &           + clund*(rop(lr,2)-ro_planlund(li,2))*amor
     &          )
      rve = roe*(   ro_inflow(li,3)
     &           + clund*(rop(lr,3)-ro_planlund(li,3))*amor
     &          )
      rwe = roe*(   ro_inflow(li,4)
     &           + clund*(rop(lr,4)-ro_planlund(li,4))*amor
     &          )

      ete = roe* ro_inflow(li,5)*param_real(CVINF)
     &     +0.5*( rue*rue + rve*rve + rwe*rwe) *roe_inv

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

