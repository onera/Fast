      lij  =       inddm( ind_loop(1) , j, k) -1
      ltij = lij - indmtr(ind_loop(1) , j, k) +1
CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
#ifdef _OPENMP4
CCCCCDIR$ VECTOR TEMPORAL
!$OMP simd aligned(dvardc,rop: CACHELINE)
#else
!DIR$ IVDEP
#endif
      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

        lt = l  - ltij
        lvo= lt
