
!rayon spectral      
        s_1 =1.0/max(sk,1.E-30)
        nx = tcx*s_1
        ny = tcy*s_1
        nz = tcz*s_1

        r_1=1./(1.+sqrt(qm1/qp1))
        

        u=qp2*r_1+qm2*(1.-r_1)
        v=qp3*r_1+qm3*(1.-r_1)
        w=qp4*r_1+qm4*(1.-r_1)
          

        f1 = abs(u*tcx+ v*tcy+ w*tcz) ! Ab ??
        f2= sqrt((tcx*tcx+tcy*tcy+tcz*tcz)*(gam1*rgp*qp5*r_1
     &  +gam1*rgp*qm5*(1.-r_1)))
        f3=0.5*wig(l+wigd)
        f4=1.-wig(l+wigd)


      !centre - tdu (ausm) - dissidaption (spectral)
        flu1 =f4*us*(r1+r2) + f3*(
     &  qn1*r1+qn2*r2) -(f4*tdu+f3*(f1+f2))
     & *(r2-r1)
        flu2 =f4*us*(rou1+rou2) + f3*(
     &  qn1*rou1+qn2*rou2) -(f4*tdu+f3*(f1+f2))
     & *(rou2-rou1)+tcx*p1p2
        flu3 =f4*us*(rov1+rov2) + f3*(
     &  qn1*rov1+qn2*rov2) -(f4*tdu+f3*(f1+f2))
     & *(rov2-rov1)+tcy*p1p2
        flu4 =f4*us*(row1+row2) + f3*(
     &  qn1*row1+qn2*row2) -(f4*tdu+f3*(f1+f2))
     & *(row2-row1)+tcz*p1p2
        flu5 =f4*us*(h1+h2) + f3*(
     &  qn1*h1+qn2*h2) -(f4*tdu+f3*(f1+f2))
     & *(h2-h1)


