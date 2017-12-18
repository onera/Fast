        r1  =qp1
        rou1=r1*qp2
        rov1=r1*qp3
        row1=r1*qp4
        p1  =r1*qp5*rgp
        h1  =gam1*p1 + .5*(rou1*qp2+rov1*qp3+row1*qp4)

        !determination etat droite: ro, roui, roe+p
        r2  =qm1
        rou2=r2*qm2
        rov2=r2*qm3
        row2=r2*qm4
        p2  =r2*qm5*rgp
        h2  =gam1*p2 + .5*(rou2*qm2+rov2*qm3+row2*qm4)
