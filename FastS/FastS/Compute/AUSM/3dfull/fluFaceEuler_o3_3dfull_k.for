c.....Metrique
 
#include  "FastS/Compute/Normale/normale_3dfull_k.for"

        nm  = l -  inck
        nm2 = l -2*inck
        np  = l +  inck

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
        vslp = v1
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm1 = qm
        qp1 = qp

        vslp = v2
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm2 = qm
        qp2 = qp

        vslp = v3
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm3 = qm
        qp3 = qp

        vslp = v4                    
#include  "FastS/Compute/Slope/o3_slope_var.for"   
        qm4 = qm                     
        qp4 = qp                     

        vslp = v5
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm5 = qm
        qp5 = qp

!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_GD.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_3dfull_k.for"

       ! modification de vitesse normale par ajout
        ! de stabilisation de type Rhie-Chow
        c   = rgp*gam*rop(l +v5)  !c^2
        son =sqrt(qn1*qn1 / c)
        tam =c3*son+sk
        tam1=max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
        u   =0.25*(qn1+qn2)-tam1*(p2-p1)
        tdu = max(abs(u),c1*sk)

        !Calcul du flux total
        p1p2= (p1+p2)*0.5


#include "FastS/Compute/Vit_ent/fluvector_3dfull_k.for"
