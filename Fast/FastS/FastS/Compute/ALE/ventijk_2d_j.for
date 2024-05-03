            lx1 = l111 + inck               ! x(i  , j  , k+1)
            lx2 = l111 + inci               ! x(i+1, j  , k  )
            lx3 = lx1  + inci               ! x(i+1, j  , k+1)

#include "FastS/Compute/ALE/target_pt.for"

            cax =  rot(1,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenz - param_real(ROT_CENTER+2))

            ventj(l ) = vtrans(1) - rot(4,3)*cay
            ventj(l2) = vtrans(2) + rot(4,3)*cax
