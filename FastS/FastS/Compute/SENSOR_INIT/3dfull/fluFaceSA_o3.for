      vslp = v6
#include  "FastS/Compute/Slope/o3_slope_var.for"
      qm6  = qm*qm1
      qp6  = qp*qp1

      flu6 = u*(qm6+qp6) -tdu*(qm6 - qp6)
