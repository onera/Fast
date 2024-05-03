      if(nitrun.le.1.and.nitcfg.eq.1) then
!$OMP BARRIER
      if(ithread.eq.param_int( IO_THREAD ) ) then

        do j=1,4
          totf(j) =0
          do i=1,Nbre_thread_actif
             totf(j) = totf(j) + tot(j,i)
          enddo
          glob(j)=glob(j)+totf(j)
#if CHECK_SPLIT > 1
          write(*,'(a,13i9)')'tot=',tot(j,:),totf(j)
#endif
          
        enddo
      endif
!$OMP BARRIER
      do j=1,4
        tot(j,ithread)=0
      enddo
      endif
