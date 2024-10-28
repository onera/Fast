        l1   = inddm( i, j, kr)
        li   = indbci(i, j) 
#include  "FastS/BC/linearization_state.for"
C ...  Implementation elsA 
        roc0 = SQRT(ro0*gam*p0)
        
        dqn  = (pext(li)-pi)/roc0
        
        rog = roi + ro0**2.*dqn/roc0

        ug = ui - dqn*tnx
        vg = vi - dqn*tny
        wg = wi - dqn*tnz
        
        rop(l,1) = rog
        rop(l,2) = ug
        rop(l,3) = vg
        rop(l,4) = wg
        rop(l,5) = pext(li)/(rog*rgp)
