       if(nitrun.le.1.and.nitcfg.eq.1) then
       inc1        = (ind_sdm(2)-ind_sdm(1)+1)*(ind_sdm(4)-ind_sdm(3)+1)
     &              *(ind_sdm(6)-ind_sdm(5)+1)
       inc2        = (ind_mjr(2)-ind_mjr(1)+1)*(ind_mjr(4)-ind_mjr(3)+1)
     &              *(ind_mjr(6)-ind_mjr(5)+1)

       inc11= inc11+inc1
       inc22= inc22+inc2

       tot(1,ithread)= tot(1,ithread)+ inc1
       tot(2,ithread)= tot(2,ithread)+ inc2
       tot(3,ithread)= tot(3,ithread)+ (ind_rhs(2)-ind_rhs(1)+1)
     &                                *(ind_rhs(4)-ind_rhs(3)+1)
     &                                *(ind_rhs(6)-ind_rhs(5)+1)
       tot(4,ithread)= tot(4,ithread)+ (ind_grad(2)-ind_grad(1)+1)
     &                                *(ind_grad(4)-ind_grad(3)+1)
     &                                *(ind_grad(6)-ind_grad(5)+1)

#if CHECK_SPLIT > 1
       if(ithread.eq.param_int( IO_THREAD).and.nitcfg.le.1) then
        write(*,*)'tot_sdm, mjr',inc1,inc11,inc2,inc22
       endif
#endif 
      endif
