
C****************************************************************************************
c***********                        DEFINITION DE W                      ****************
C****************************************************************************************

       Cpm = (gam*rgp)/(gam-1) ! Capacité calorifique massique à pression constante
       Cvm = rgp/(gam-1) ! Capacité calorifique massique à volume constant

c------------- En np
       W1_np = rop(np+v1)
       W2_np = rop(np+v1) * rop(np+v2)
       W3_np = rop(np+v1) * rop(np+v3)
       rho_e_cin_np = 0.5*(W2_np**2 + W3_np**2)
       W5_np = Cvm * W1_np * rop(np+v5) + rho_e_cin_np/W1_np

c------------- En l
       W1_l = rop(l+v1)
       W2_l = rop(l+v1) * rop(l+v2)
       W3_l = rop(l+v1) * rop(l+v3)
       rho_e_cin_l = 0.5*(W2_l**2 + W3_l**2)
       W5_l = Cvm * W1_l * rop(l+v5) + rho_e_cin_l/W1_l


C****************************************************************************************
c***********                      CALCUL de dpP_dpW                      ****************
C****************************************************************************************

c------------- En np
       dpRop1_dpW1_np = 1.

       dpRop2_dpW1_np = -W2_np/(W1_np**2)
       dpRop2_dpW2_np = 1./W1_np

       dpRop3_dpW1_np = -W3_np/(W1_np**2)
       dpRop3_dpW3_np = 1./W1_np

       dpRop5_dpW1_np = 1./Cvm * (- W5_np/(W1_np**2) +
     &                  2 * rho_e_cin_np/(W1_np**3))
       dpRop5_dpW2_np = - 1./Cvm * W2_np/(W1_np**2)
       dpRop5_dpW3_np = - 1./Cvm * W3_np/(W1_np**2)
       dpRop5_dpW5_np = 1./ (Cvm * W1_np)

c------------- En l
       dpRop1_dpW1_l = 1.

       dpRop2_dpW1_l = -W2_l/(W1_l**2)
       dpRop2_dpW2_l = 1./W1_l

       dpRop3_dpW1_l = -W3_l/(W1_l**2)
       dpRop3_dpW3_l = 1./W1_l

       dpRop5_dpW1_l = 1./Cvm * (- W5_l/(W1_l**2) +
     &                 2 * rho_e_cin_l/(W1_l**3))
       dpRop5_dpW2_l = - 1./Cvm * W2_l/(W1_l**2)
       dpRop5_dpW3_l = - 1./Cvm * W3_l/(W1_l**2)
       dpRop5_dpW5_l = 1./(Cvm * W1_l)

