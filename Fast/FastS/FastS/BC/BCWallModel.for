       ventx= ventijk(ldp,    1  )*c_ale*mobile_coef          ! rotational/translational velocity
       venty= ventijk(ldp,    2  )*c_ale*mobile_coef          ! rotational/translational velocity 
       ventz= ventijk(ldp,kc_vent)*ck_vent*c_ale*mobile_coef  ! rotational/translational velocity

       !etat aero externe CL
       text = rop(loo , 5) 
       uext = sqrt( (rop(loo, 2)-ventx)**2
     &             +(rop(loo, 3)-venty)**2
     &             +(rop(loo, 4)-ventz)**2 )
       uext=max(uext,1e-15)
       tt_ext= text + 0.5*uext**2*cpinv

       unorm= sqrt( (rop(ldjr, 2)-ventx)**2
     &             +(rop(ldjr, 3)-venty)**2
     &             +(rop(ldjr, 4)-ventz)**2 )

       rop(l , 2)=c7*(rop(ldl , 2) -rop(ldnm , 2) ) +rop(ldnp , 2)
       rop(l , 3)=c7*(rop(ldl , 3) -rop(ldnm , 3) ) +rop(ldnp , 3)
       rop(l , 4)=c7*(rop(ldl , 4) -rop(ldnm , 4) ) +rop(ldnp , 4)


       !0.70: coeff optimisable lier a bflwallmodel
        twall = text + 0.5*0.70*uext**2*cpinv
c      twall = text + 0.5*0.90*(uext**2-unorm**2)*cpinv


       rop(l , 5) = 2* twall - rop(ldnm , 5)

       !cutoff robustesse
       rop(l , 5) = min(rop(l , 5), tt_ext)
       rop(l , 5) = max(rop(l , 5), text)


       !dpdn=0 (a partir point image)
       rop(l , 1) = rop(m, 1)*rop(m, 5)/rop(l , 5)

c       if(( j.le.2.or.j.ge.155).and.i.ge.155) then
c           write(*,'(a,8f8.2,2i4)')'bvbs',tt_ext, text, 
c     &   twall, rop(ldl+v5), rop(ldnm +v5),rop(l+v5),  unorm, uext,i,j
c       endif
c            if(i.eq.100.or.i.eq.120) then
c           write(*,'(a,8f8.2,i3)')'bvbs',tt_ext, text, 
c     &    twall, rop(ldl+v5), rop(ldnm +v5),rop(l+v5),  unorm, uext,i
c            endif
