#ifndef E_SCALAR_COMPUTER
             do k = ind_loop(5), ind_loop(6)
             do i = ind_loop(1), ind_loop(2)

               l = inddm( i,j,k)
#else
             do k = ind_loop(5), ind_loop(6)
               lij  =       inddm( ind_loop(1) , j, k)
               ltij = lij - indmtr(ind_loop(1) , j, k)
#ifdef _OPENMP4
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
                 do l = lij, lij +  ind_loop(2)- ind_loop(1) 

#endif 
