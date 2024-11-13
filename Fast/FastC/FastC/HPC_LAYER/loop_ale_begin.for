#ifndef E_SCALAR_COMPUTER
#include   "FastS/Compute/loopGpu_begin.for"
            lt =  indmtr(i, j, k)
            lv =  indven(i, j, k)
            lx =  indcg( i, j, k)
#else
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)
                ltij = lij - indmtr(ind_loop(1) , j, k)
                lvij = lij - indven(ind_loop(1) , j, k)
                lxij = lij - indcg( ind_loop(1) , j, k)
#ifdef _OPENMP4
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
                do l = lij, lij +  ind_loop(2)- ind_loop(1)
                    lt = l  - ltij
                    lv = l  - lvij
                    lx = l  - lxij
#endif 
                    lvo= lt
