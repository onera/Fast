#ifndef E_SCALAR_COMPUTER
             do j = ind_loop(3), ind_loop(4)
             do i = ind_loop(1), ind_loop(2)

               l = inddm( i,j,k)
               lt= indmtr(i,j,k)
#else
             do j = ind_loop(3), ind_loop(4)
               lij  =       inddm( ind_loop(1) , j, k)
               ltij = lij - indmtr(ind_loop(1) , j, k)
#ifdef _OPENMP4
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
                 do l = lij, lij +  ind_loop(2)- ind_loop(1) 

               lt= l-ltij
#endif 
               lvo=lt 
