#include  "Fast/BC/BCFarfield_firstrank.for"

      if(eigv1.lt.0.) then
            rop(ir,6) = nue
      else
            rop(ir,6) = rop(il,6)
      endif
