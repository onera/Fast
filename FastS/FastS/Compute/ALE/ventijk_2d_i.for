            lx1 = l111 + incj
            lx2 = l111 + inck
            lx3 = lx1  + inck

#include "FastS/Compute/ALE/target_pt.for"

            cax =  rot(1,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenz - param_real(ROT_CENTER+2))

            venti(l ) = vtrans(1) - rot(4,3)*cay
            venti(l2) = vtrans(2) + rot(4,3)*cax 

            !write(*,*)venti(l ),venti(l2 ),j,l-lij
