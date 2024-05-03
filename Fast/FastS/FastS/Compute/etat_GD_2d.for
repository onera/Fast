        r1  =qp1        ! rho
        rou1=r1*qp2     ! rho*u
        rov1=r1*qp3     ! rho*v
        p1  =r1*qp5*rgp ! rho*T*R
        h1  =gam1*p1 + .5*(rou1*qp2+rov1*qp3) !total enthalpy

        !determination etat droite: ro, roui, roe+p
        r2  =qm1
        rou2=r2*qm2
        rov2=r2*qm3
        p2  =r2*qm5*rgp
        h2  =gam1*p2 + .5*(rou2*qm2+rov2*qm3)
