
C****************************************************************************************
c***********                        DEFINITION DE W                      ****************
C****************************************************************************************

       Cpm = (gam*rgp)/(gam-1) ! Capacité calorifique massique à pression constante
       Cvm = rgp/(gam-1) ! Capacité calorifique massique à volume constant

c------------- En np
       W1_nm2 = rop(nm2+v1)
       W2_nm2 = rop(nm2+v1) * rop(nm2+v2)
       W3_nm2 = rop(nm2+v1) * rop(nm2+v3)
       W4_nm2 = rop(nm2+v1) * rop(nm2+v4)
       rho_e_cin_nm2 = 0.5*(W2_nm2**2 + W3_nm2**2 + W4_nm2**2)
       W5_nm2 = Cvm * W1_nm2 * rop(nm2+v5) + rho_e_cin_nm2/W1_nm2

c------------- En l
       W1_nm = rop(nm+v1)
       W2_nm = rop(nm+v1) * rop(nm+v2)
       W3_nm = rop(nm+v1) * rop(nm+v3)
       W4_nm = rop(nm+v1) * rop(nm+v4)
       rho_e_cin_nm = 0.5*(W2_nm**2 + W3_nm**2 + W4_nm**2)
       W5_nm = Cvm * W1_nm * rop(nm+v5) + rho_e_cin_nm/W1_nm


C****************************************************************************************
c***********                      CALCUL de dpP_dpW                      ****************
C****************************************************************************************

c------------- En np
       dpRop1_dpW1_nm2 = 1.

       dpRop2_dpW1_nm2 = -W2_nm2/(W1_nm2**2)
       dpRop2_dpW2_nm2 = 1./W1_nm2

       dpRop3_dpW1_nm2 = -W3_nm2/(W1_nm2**2)
       dpRop3_dpW3_nm2 = 1./W1_nm2

       dpRop4_dpW1_nm2 = -W4_nm2/(W1_nm2**2)
       dpRop4_dpW4_nm2 = 1./W1_nm2

       dpRop5_dpW1_nm2 = 1./Cvm * (- W5_nm2/(W1_nm2**2) +
     &                   2 * rho_e_cin_nm2/(W1_nm2**3))
       dpRop5_dpW2_nm2 = - 1./Cvm * W2_nm2/(W1_nm2**2)
       dpRop5_dpW3_nm2 = - 1./Cvm * W3_nm2/(W1_nm2**2)
       dpRop5_dpW4_nm2 = - 1./Cvm * W4_nm2/(W1_nm2**2)
       dpRop5_dpW5_nm2 = 1./ (Cvm * W1_nm2)

c------------- En l
       dpRop1_dpW1_nm = 1.

       dpRop2_dpW1_nm = -W2_nm/(W1_nm**2)
       dpRop2_dpW2_nm = 1./W1_nm

       dpRop3_dpW1_nm = -W3_nm/(W1_nm**2)
       dpRop3_dpW3_nm = 1./W1_nm

       dpRop4_dpW1_nm = -W4_nm/(W1_nm**2)
       dpRop4_dpW4_nm = 1./W1_nm

       dpRop5_dpW1_nm = 1./Cvm * (- W5_nm/(W1_nm**2) +
     &                  2 * rho_e_cin_nm/(W1_nm**3))
       dpRop5_dpW2_nm = - 1./Cvm * W2_nm/(W1_nm**2)
       dpRop5_dpW3_nm = - 1./Cvm * W3_nm/(W1_nm**2)
       dpRop5_dpW4_nm = - 1./Cvm * W4_nm/(W1_nm**2)
       dpRop5_dpW5_nm = 1./(Cvm * W1_nm)

