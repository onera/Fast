      do k = ind_loop(5), ind_loop(6)
      do j = ind_loop(3), ind_loop(4)
        lij  = inddm(ind_loop(1) , j, k)-1
#ifdef _OPENMP4
CCCcDEC$ PREFETCH rop_1
CCC  !$OMP SIMD NONTEMPORAL(rop_n1)
CC!DIR$ ASSUME_ALIGNED drodm: 32
CC!DIR$ ASSUME_ALIGNED coe: 32
CC!DIR$ ASSUME (mod(param_int(NDIMDX), VECLENGTH) .eq. 0)
CC!DIR$ ASSUME (mod(lij,4) .eq. 0)
!$OMP SIMD
#else
CDIR$ IVDEP
#endif
        do l = lij+1, lij+1 +  ind_loop(2) - ind_loop(1)
