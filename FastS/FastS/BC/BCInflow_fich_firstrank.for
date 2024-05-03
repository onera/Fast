      roe     = dens(li)
      roe_inv = 1./roe       

      rue=ux(li)*dens(li)
      rve=uy(li)*dens(li)
      rwe=uz(li)*dens(li)

      ete = dens(li)* temp(li)*param_real(CVINF)
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
      !modif ivan pour supprimer onde parasite acoustique
      !rop(l,1) =  roe
      !rop(l,2) =  ux(li) 
      !rop(l,3) =  uy(li)
      !rop(l,4) =  uz(li) 
      !rop(l,5) =  rop(l1,5)*rop(l1,1)*roe_inv 
