        ! Etat (qm) a l'interface droite 

        qm1=  c4*rop(l +v1)  +  c5*rop(nm +v1)  +  c6*rop(np +v1)
        qm2=  c4*rop(l +v2)  +  c5*rop(nm +v2)  +  c6*rop(np +v2)
        qm3=  c4*rop(l +v3)  +  c5*rop(nm +v3)  +  c6*rop(np +v3)
        qm5=  c4*rop(l +v5)  +  c5*rop(nm +v5)  +  c6*rop(np +v5)

        ! Etat (qp) a l'interface gauche

        qp1=  c4*rop(nm +v1)  +  c6*rop(nm2 +v1) +  c5*rop(l +v1)
        qp2=  c4*rop(nm +v2)  +  c6*rop(nm2 +v2) +  c5*rop(l +v2)
        qp3=  c4*rop(nm +v3)  +  c6*rop(nm2 +v3) +  c5*rop(l +v3)
        qp5=  c4*rop(nm +v5)  +  c6*rop(nm2 +v5) +  c5*rop(l +v5)

