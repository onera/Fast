      !Pentes annulees pour toute temperature inferieure a la temperature Tmin
      roff=.5*max(0.,min(rop(l +v5)*Tmin_1 - .5,1.))

      ! Pente densite
      delm=rop(l  +v1) - rop(nm +v1)
      delp=rop(np +v1) - rop(l  +v1)
      delq=rop(nm +v1) - rop(nm2+v1)


      slp =roff*avmin(delm,delp)
      slq =roff*avmin(delm,delq)

      qp1   =rop(nm +v1) +slq
      qm1   =rop(l  +v1) -slp

      ! Pente velocityX
      delm=rop(l  +v2) - rop(nm +v2)
      delp=rop(np +v2) - rop(l  +v2)
      delq=rop(nm +v2) - rop(nm2+v2)


      slp =roff*avmin(delm,delp)
      slq =roff*avmin(delm,delq)

      qp2  =rop(nm +v2) +slq
      qm2  =rop(l  +v2) -slp

      ! Pente velocityY
      delm=rop(l  +v3) - rop(nm +v3)
      delp=rop(np +v3) - rop(l  +v3)
      delq=rop(nm +v3) - rop(nm2+v3)

      slp =roff*avmin(delm,delp)
      slq =roff*avmin(delm,delq)

      qp3  =rop(nm +v3) +slq
      qm3  =rop(l  +v3) -slp

      ! Pente temperature
      delm=rop(l  +v5) - rop(nm +v5)
      delp=rop(np +v5) - rop(l  +v5)
      delq=rop(nm +v5) - rop(nm2+v5)

      slp =roff*avmin(delm,delp)
      slq =roff*avmin(delm,delq)

      qp5  =rop(nm +v5) +slq
      qm5  =rop(l  +v5) -slp
