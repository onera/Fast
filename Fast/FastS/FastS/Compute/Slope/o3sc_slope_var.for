        ! qm: right state,  qp: left state
        qm=  c4*rop(l +vslp) +  c5*rop(nm  +vslp) +  c6*rop(np +vslp)
        qp=  c4*rop(nm+vslp) +  c6*rop(nm2 +vslp) +  c5*rop(l  +vslp)
        qm=     rop(l +vslp)*(1.-c) + qm*c 
        qp=     rop(nm+vslp)*(1.-c) + qp*c

