#ifndef E_SCALAR_COMPUTER
             do j = ind_loop(3), ind_loop(4)
             do i = ind_loop(1), ind_loop(2)
#else
             do j = ind_loop(3), ind_loop(4)
#ifdef _OPENMP4
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
              do i = ind_loop(1), ind_loop(2)
#endif 
               l = inddm( i,j,k)
