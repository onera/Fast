c.....Metrique
 
#include  "FastS/Compute/Normale/normale_3dhomo_j.for"

        nm  = l -  incj
        nm2 = l -2*incj
        np  = l +  incj

#include  "FastS/Compute/SENSOR_INIT/wiggle_j.for"

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
#include "FastS/Compute/Vit_ent/qn_ale_3dhomo_j.for"

        ! modification de vitesse normale par ajout
        ! de stabilisation de type Rhie-Chow
        u   = 0.25*(qn1+qn2)- c2*sj*(p2-p1)*wig( l+ wig_j)

        tdu = max(abs(u),c1*sj)*wig( l+ wig_j)

        !Calcul du flux total
        p1p2= (p1+p2)*0.5

#include "FastS/Compute/Vit_ent/fluvector_ale_3dhomo_j.for"
