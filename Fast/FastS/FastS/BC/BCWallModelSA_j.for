            ldjr = inddm(  i, jr        , k )
            ldl  = inddm(  i, jr +sens  , k )
            ldnp = inddm(  i, jr +sens*2, k )
            ldnm = ldjr
 
            m    = ldjr  - exchange
            loo  = ldjr  - exchange*shift_loo
           !remplissage ghost 5 variables
#include     "FastS/BC/BCWallModel.for"
           !write(*,*)'rop',rop(l,1),rop(l,2),rop(l,3),rop(l,5),j,i

            iadrf = ldjr
            mvent = ldp  - shiftvent
           !calcul de utau
#include "FastS/BC/BFWALLMODEL_2indice.for"
           !write(*,*)'utau',utau
           !forcage nutilde premiere cellule reelle: nut= kappa y+
             ltg= ldjr
#include  "FastS/BC/BCWallExchange_nutilde.for"
             !!remplissage ghost pour pente gauche = pentre droite
             rop(l , 6)=c7*(rop(ldl , 6) -rop(ldnm , 6) ) +rop(ldnp , 6)
