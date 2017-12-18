      ltij  =       indmtr( ind_loop(1) , j, k) -1
CC    !DIR$ ASSUME (mod(ltij,4) .eq. 0)
#ifdef _OPENMP4
CCCC!$OMP simd aligned(ti,tj,tk,ti0,tj0,tk0: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
      do lt = ltij+1, ltij+1 + ind_loop(2) - ind_loop(1)
