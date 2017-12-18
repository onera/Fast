      !Pentes annulees pour toute temperature inferieure a la temperature Tmin
      roff=.5*max(0.,min(rop(l +v5)*Tmin_1 - .5,1.))

      ! qm: right state,  qp: left state
      delm=rop(l  +vslp) - rop(nm +vslp)  !i+1 - i 
      delp=rop(np +vslp) - rop(l  +vslp)  !i+2 - i+1
      delq=rop(nm +vslp) - rop(nm2+vslp)  !i   - i-1


      slp =roff*avmin(delm, delp)
      slq =roff*avmin(delq, delm)

      !right state
      qm   =rop(l  +vslp) -slp

      !left state
      qp   =rop(nm +vslp) +slq
