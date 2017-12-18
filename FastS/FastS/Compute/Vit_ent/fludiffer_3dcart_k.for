        du =(qm2-qp2)*sk
        dv =(qm3-qp3)*sk
        dw =(qm4-qp4)*sk
        dp = p2-p1
        dqn=qn2-qn1  !Sdqn

c.....Evaluation des valeurs propres
        s_1 =1./max(sk,1.e-30)
        nz = tcz*s_1

        qn = w*nz  

        f1 = qn
        f4 = psiroe*(abs(f1) + c)  !delta

        f2 =.5+sign(.5,abs(f1)-f4)
        f3 =(f2-1.)*(f1*f1/f4 + f4)*0.5  -f2*abs(f1)
 
        f1  = (r2-r1 - dp/(c*c))*sk !  a123
 
        flu1= f3*( f1                                   )
        flu2= f3*( f1 * u + r * (du         )           )
        flu3= f3*( f1 * v + r * (dv         )           )
        flu4= f3*( f1 * w + r * (dw - nz*dqn)           )
        flu5= f3*( f1 * q + r * (u*du+v*dv+w*dw-qn*dqn) )

        f1 = qn+c
        f2 =.5+sign(.5,abs(f1)-f4)
        f3 =(f2-1.)*(f1*f1/f4 + f4)*0.5  -f2*abs(f1)
 
        f1 = (dp+r*c*dqn*s_1) / (2.*c*c)*sk*f3  !lambda4*S*a4

        flu1 = flu1 + f1
        flu2 = flu2 + f1 * (u        )
        flu3 = flu3 + f1 * (v        )
        flu4 = flu4 + f1 * (w + nz* c)
        flu5 = flu5 + f1 * (h + qn *c)

        f1 = qn-c
        f2 =.5+sign(.5,abs(f1)-f4)
        f3 =(f2-1.)*(f1*f1/f4 + f4)*0.5  -f2*abs(f1)

        f1 = (dp-r*c*dqn*s_1) / (2.*c*c)*f3*sk !lambda5*S*a5
 
        flu1 = flu1 + f1 
        flu2 = flu2 + f1 * (u        )
        flu3 = flu3 + f1 * (v        )
        flu4 = flu4 + f1 * (w - nz* c)
        flu5 = flu5 + f1 * (h - qn* c)

        rou1=r1*qn1
        rou2=r2*qn2
        p1p2=p1+p2
 
        flu1 = (flu1 + rou1     + rou2                 )*0.5
        flu2 = (flu2 + rou1*qp2 + rou2*qm2             )*0.5
        flu3 = (flu3 + rou1*qp3 + rou2*qm3             )*0.5
        flu4 = (flu4 + rou1*qp4 + rou2*qm4 + tcz * p1p2)*0.5
        flu5 = (flu5 + rou1*h1  + rou2*h2              )*0.5
