            ldjr = inddm(  ir        , j, k )
            ldl  = inddm(  ir +sens_int  , j,  k )
            ldnp = inddm(  ir +sens_int*2, j,  k )
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
