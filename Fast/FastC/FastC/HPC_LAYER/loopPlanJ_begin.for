#ifndef E_SCALAR_COMPUTER
             do k = ind_loop(5), ind_loop(6)
             do i = ind_loop(1), ind_loop(2)

#else
             do k = ind_loop(5), ind_loop(6)
#ifdef _OPENMP4
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
             do i = ind_loop(1), ind_loop(2)
#endif 

               l = inddm( i,j,k)
