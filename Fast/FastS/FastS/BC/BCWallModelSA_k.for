            ldjr = inddm( i, j,  kr         )
            ldl  = inddm( i, j,  kr +sens_int   )
            ldnp = inddm( i, j,  kr +sens_int*2 )
            ldnm = ldjr
 
            m    = ldjr  - exchange
            loo  = ldjr  - exchange*shift_loo
           !remplissage ghost 5 variables
#include     "FastS/BC/BCWallModel.for"

            iadrf = ldjr
            mvent = ldp  - shiftvent
           !calcul de utau
#include "FastS/BC/BFWALLMODEL_2indice.for"
           !forcage nutilde premiere cellule reelle: nut= kappa y+
             ltg= ldjr
#include  "FastS/BC/BCWallExchange_nutilde.for"
             !!remplissage ghost pour pente gauche = pentre droite
             rop(l , 6)=c7*(rop(ldl , 6) -rop(ldnm , 6) ) +rop(ldnp , 6)
