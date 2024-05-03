c.....Metrique
 
#include  "FastS/Compute/Normale/normale_2d_j.for"

        nm  = l -  incj
        nm2 = l -2*incj
        nm3 = l -3*incj
        np  = l +  incj
        np2 = l +2*incj

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
        c    = wig(l+sl_j)     !SC only
        vslp = v1
#include  "FastS/Compute/Slope/o5sc_slope_var.for"
        qm1 = qm
        qp1 = qp

        vslp = v2
#include  "FastS/Compute/Slope/o5sc_slope_var.for"
        qm2 = qm
        qp2 = qp

        vslp = v3
#include  "FastS/Compute/Slope/o5sc_slope_var.for"
        qm3 = qm
        qp3 = qp


        vslp = v5
#include  "FastS/Compute/Slope/o5sc_slope_var.for"
        qm5 = qm
        qp5 = qp

!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_GD_2d.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_2d_j.for"

       ! modification de vitesse normale par ajout
        ! de stabilisation de type Rhie-Chow
        c   = rgp*gam*rop(l +v5)  !c^2
        son =sqrt(qn1*qn1 / c)
        tam =c3*son+sj
        tam1=max(0.,tam)*c2 ! fct amortissement: c3*Mach+1
        u   =0.25*(qn1+qn2)-tam1*(p2-p1)
        tdu = max(abs(u),c1*sj)

        !Calcul du flux total
        p1p2= (p1+p2)*0.5


#include "FastS/Compute/Vit_ent/fluvector_2d_j.for"
