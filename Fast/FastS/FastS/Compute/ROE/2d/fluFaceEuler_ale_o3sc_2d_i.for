c.....Metrique
 
#include  "FastS/Compute/Normale/normale_2d_i.for"

        nm  = l -  inci
        nm2 = l -2*inci
        np  = l +  inci

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
        c    = wig(l+sl_i)     !SC only
        vslp = v1
#include  "FastS/Compute/Slope/o3sc_slope_var.for"
        qm1 = qm
        qp1 = qp

        vslp = v2
#include  "FastS/Compute/Slope/o3sc_slope_var.for"
        qm2 = qm
        qp2 = qp

        vslp = v3
#include  "FastS/Compute/Slope/o3sc_slope_var.for"
        qm3 = qm
        qp3 = qp


        vslp = v5
#include  "FastS/Compute/Slope/o3sc_slope_var.for"
        qm5 = qm
        qp5 = qp

!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_roe_GD_2d.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_ale_2d_i.for"

#include "FastS/Compute/Vit_ent/fludiffer_ale_2d_i.for"
