      lij  =       inddm( ind_loop(1) , j, k) -1
      ltij = lij - indmtr(ind_loop(1) , j, k) +1
CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
#ifdef _OPENMP4
CCCCCDIR$ VECTOR TEMPORAL
CCCC!$OMP simd aligned(drodm,rop,xmut: CACHELINE)
CC!DIR$ ASSUME_ALIGNED drodm: 32
CC!DIR$ ASSUME (mod(lij,4).eq.0)
CC!DIR$ ASSUME (mod(inck,   VECLENGTH).eq.0)
CC!DIR$ ASSUME (mod(incj,   VECLENGTH).eq.0)
CC!DIR$ ASSUME (mod(param_int(NDIMDX), VECLENGTH) .eq. 0)
!$OMP SIMD 
#else
!DIR$ IVDEP
#endif
CCCC!DIR$ DISTRIBUTE POINT
      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

        lt = l  - ltij
        lvo= lt
