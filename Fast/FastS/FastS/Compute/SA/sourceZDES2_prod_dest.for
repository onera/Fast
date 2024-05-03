!!  amulam :: μ kinematic viscosity based on Sutherland
#include  "FastS/Compute/mulam.for"
!! anutild :: ν turbulent SA nu tilde  || chi:: ν/v(kinematic viscosity)
#include  "FastS/Compute/SA/chi.for"
          fvv1 = fv1(chi)
#include  "FastS/Compute/SA/xmut.for"
!!ad1 :: d_w wall distance        
          ad1      = max(dlng(l) , 1.e-27)
          
          ra       =(xmut(l)/rop(l ,1))
     1                  /(SA_CKARM*SA_CKARM*ad1*ad1*auijuij**0.5)  ! r_d = (\{nu}_t+\{nu})/(\sqrt{U_{i,j}^2}*\kappa^2*d_w^2)
          variable2=512.*ra*ra*ra                                      ! (8*r_d)^3
          variable2=min(variable2,15.)                                 ! cut tanh to 1
          fa       =1. - tanh(variable2)                               ! f_d = 1-tanh((8*r_d)^3)
          testfa   =0.5+ sign(.5, 0.8-fa)                              ! if fd>fd0: 0.5+-0.5=0 || if fd<fd0: 0.5+0.5=1
          adelta1  =testfa*sph2 +(1.-testfa)*adelta1                   ! echelle caractéristique maillage

          !choix distance
          adtild1  =ad1-fa*max(ad1-SA_RCTEDES*adelta1,0.)              !d_{DES}^{II}=d_w - f_d*max(d_w-C_{DES}\Delta_{DES}^{II},0)
          dist     =max(adtild1,1.e-27)                                ! make sure we do not divide by 0

          !CALCUL DU TERME DE PRODUCTION
          f2       = fv2(chi,fvv1)
#include  "FastS/Compute/SA/zdes2_zdes3_common.for"
          fwg       = fw(g)
#include  "FastS/Compute/SA/zdes2_zdes3_common_prt2.for"

