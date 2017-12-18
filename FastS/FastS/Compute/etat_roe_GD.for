        r1=qp1
        p1=r1*qp5*rgp
        h1=gam1*qp5*rgp + .5*(qp2*qp2+qp3*qp3+qp4*qp4)

        !determination etat droite: ro, roui, roe+p
        r2 =qm1
        p2=r2*qm5*rgp
        h2=gam1*qm5*rgp + .5*(qm2*qm2+qm3*qm3+qm4*qm4)

!determination etat moyenne roe
        qp1=sqrt(r1)
        qm1=sqrt(r2)
        r  = qp1*qm1
        r_1= 1./(qp1+qm1)
        u =(qp1*qp2+qm1*qm2)*r_1
        v =(qp1*qp3+qm1*qm3)*r_1
        w =(qp1*qp4+qm1*qm4)*r_1
        h =(qp1*h1 +qm1*h2 )*r_1
        q =.5*(u*u+v*v+w*w)
        c =sqrt((gam-1.)*abs(h-q))
