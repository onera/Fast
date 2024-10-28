#ifndef E_SCALAR_COMPUTER
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)

#else
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
               !lij  =       inddm( ind_loop(1) , j, k)
               !ltij = lij - indmtr(ind_loop(1) , j, k)
C#ifdef _OPENMP4
C!$OMP SIMD 
C#else
C!DIR$ IVDEP
C#endif

#endif 
               l = inddm(i,j,k)
