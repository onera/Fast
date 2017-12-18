            lx1 = l111 + incj
            lx2 = l111 + inck
            lx3 = lx1  + inck

            cenx = .25*( x(l111) + x(lx2) + x(lx3) + x(lx1) )
            ceny = .25*( y(l111) + y(lx2) + y(lx3) + y(lx1) )
            cenz = .25*( z(l111) + z(lx2) + z(lx3) + z(lx1) )

            cax =  rot(1,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenz - param_real(ROT_CENTER+2))

            caz =  rot(3,1)*(cenx - param_real(ROT_CENTER  ))
     &           + rot(3,2)*(ceny - param_real(ROT_CENTER+1))
     &           + rot(3,3)*(cenz - param_real(ROT_CENTER+2))

            venti(l ) = vtrans(1) + rot(4,2)*caz - rot(4,3)*cay
            venti(l2) = vtrans(2) + rot(4,3)*cax - rot(4,1)*caz
            venti(l3) = vtrans(3) + rot(4,1)*cay - rot(4,2)*cax
