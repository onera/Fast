      do k = ind_loop(5), ind_loop(6)
      do j = ind_loop(3), ind_loop(4)
        lij  = inddm(ind_loop(1) , j, k)
#ifdef _OPENMP4
cDEC$ PREFETCH rop_1
CDIR$ VECTOR NONTEMPORAL (rop_1)
CC!$OMP simd aligned(drodm,rop,rop_1,coe: CACHELINE)
!$OMP simd
#else
CDIR$ IVDEP
#endif
        do l = lij, lij +  ind_loop(2) - ind_loop(1)
