      lij  =       indven( ind_loop(1) , j, k) -1
      lxij = lij - indcg(  ind_loop(1) , j, k) +1      
CC    !DIR$ ASSUME (mod(lij,4) .eq. 0)
#ifdef _OPENMP4
CCCC!$OMP simd aligned(venti,ventj,ventk: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
      do l = lij+1, lij+1 + ind_loop(2) - ind_loop(1)

            l2  = l +  v2ven
            l3  = l +  v3ven
            l111= l  - lxij
