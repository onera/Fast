           ! minmod nutilde 
      vslp = v6
#include  "Slope/o3_slope_var.for"
      qm6  = qm*r2*qn2
      qp6  = qp*r1*qn1

      flu6 = 0.5*(qm6 + qp6 +sign(1.,flu1)*(qp6-qm6))
