      l1 =inddm(i,jr,k)
#include  "FastS/BC/BCFarfield.for"
      if(eigv1.lt.0.) then
            rop(l,6) = nue
      else
            rop(l,6) = rop(l1,6)
      endif

