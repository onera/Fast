      lij  =       inddm( ind_loop(1) , j, k) -1
#ifdef _OPENMP4
!$OMP simd 
#else
!DIR$ IVDEP
#endif
      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)
