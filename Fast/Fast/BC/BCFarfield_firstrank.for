#include  "Fast/BC/non_reflection.for"

      roinv = 1./svar1

      u  = u + r*svar2*roinv
      v  = v + r*svar3*roinv
      w  = w + r*svar4*roinv


      rop(ir,1) =  svar1
      rop(ir,2) =  u 
      rop(ir,3) =  v 
      rop(ir,4) =  w 
      rop(ir,5) =  (qvar5*roinv - .5*(u*u+v*v+w*w))*cvinv

