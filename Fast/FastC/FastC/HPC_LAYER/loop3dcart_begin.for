#ifndef E_SCALAR_COMPUTER
#include   "FastS/Compute/loop3dcartGpu_begin.for"
#else
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)
#ifdef _OPENMP4
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
                do l = lij, lij +  ind_loop(2)- ind_loop(1)
#endif 
