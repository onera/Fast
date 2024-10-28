        li   = indbci(j,  k )
        l1   = l + inci

        !Init Newton
#include "FastS/BC/BCInflow_newton_1_firstrank.for"

        !resolution Newton             
        residug  = 1.e+20
        nitnwt   = 0
        DO WHILE (residug .GT. tolnewton .AND.nitnwt .LT. newtonmax)
                 
            nitnwt  = nitnwt+1
            residug = 0.
#include    "FastS/BC/BCInflow_newton_2_firstrank.for"
        ENDDO ! End Newton
               
         !Mise a jour first rank
#include "FastS/BC/BCInflow_newton_3_firstrank.for"
