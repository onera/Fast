      l1 =inddm(i,j,kr)
#include  "FastS/BC/BCFarfield.for"
      if(eigv1.lt.0.) then
            rop(l,6) = nue
      else
            rop(l,6) = rop(l1,6)
      endif

