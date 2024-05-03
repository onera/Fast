        ! qm: right state,  qp: left state
        qm= ( c50*rop(l +vslp) +c51*rop(nm  +vslp)+c52*rop(np+vslp)
     &  +c53*rop(nm2 +vslp) + c54*rop(np2 + vslp))


        qp=  (c50*rop(nm+vslp) +c52*rop(nm2 +vslp)+c51*rop(l +vslp) 
     &  +c53*rop(np+vslp) + c54*rop(nm3+vslp))

