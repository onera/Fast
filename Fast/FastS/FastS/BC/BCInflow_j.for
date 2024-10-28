      l1 =inddm(i, jr,k)
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
      !rop(l,5) =  rop(l1,5)*rop(l1,1)*roe_inv 
