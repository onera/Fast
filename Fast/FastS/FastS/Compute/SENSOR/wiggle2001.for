c***********************************************************************
c     $Date: 2023-07-30 09:20:20 
c     $Revision: 1 $
c     $Author: Ivan Mary $
c***********************************************************************
      subroutine wiggle2001(ndom, param_int, param_real,
     &                      ind_loop, rop, wig)
c***********************************************************************
c
c***********************************************************************
      implicit  none

      REAL_E souszero, cutoff
      parameter(souszero=-1e-12)
      parameter(cutoff=1e-10)

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), param_int(0:*)

      REAL_E rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E wig( param_int(NDIMDX)*3)
      REAL_E param_real(0:*)
    
c Var loc 
      INTEGER_E incmax,i,j,k,l,inci,incj,inck,lij,ltij,lt,lvo,np2,np3,
     & nm,nm2,np, v1,v2,v3,v4,v5,icorr,jcorr,kcorr

      REAL_E qm1,qp1,qm2,qp2,qm3,qp3,qm4,qp4,qm5,qp5,f1,f2,f3,f4,f5,test
      REAL_E c3, c5,roref2_inv,vref2_inv,tref2_inv

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

      icorr=0
      jcorr=0
      kcorr=0
      if(ind_loop(2).eq.param_int(IJKV))   icorr=1
      if(ind_loop(4).eq.param_int(IJKV+1)) jcorr=1
      if(ind_loop(6).eq.param_int(IJKV+2)) kcorr=1

      roref2_inv = 1./(param_real(ROINF)*param_real(ROINF))
      vref2_inv  = 1./(param_real(VINF)*param_real(VINF))
      tref2_inv  = 1./(param_real(TINF)*param_real(TINF))

      IF(param_int(ITYPZONE).eq.3)  THEN

#include "FastC/HPC_LAYER/loop_begin.for"                  
           nm  = l -  incj
           nm2 = l -2*incj
           np  = l +  incj
#include   "FastS/Compute/SENSOR/wiggle2d_j.for"
           nm  = l -  inci
           nm2 = l -2*inci
           np  = l +  inci
#include   "FastS/Compute/SENSOR/wiggle2d_i.for"
#include "FastC/HPC_LAYER/loop_end.for" 
                 
          if(icorr.eq.1) then !flux manquant en I
            i   = ind_loop(2) + 1
#include    "FastC/HPC_LAYER/loopPlanI_begin.for"
               nm  = l -  inci
               nm2 = l -2*inci
               np  = l +  inci
#include       "FastS/Compute/SENSOR/wiggle2d_i.for"
#include    "FastC/HPC_LAYER/loopPlan_end.for"
          endif !
         !complement jmax
         If(jcorr.eq.1) then
           j    = ind_loop(4)+1
#include    "FastC/HPC_LAYER/loopPlanJ_begin.for"
             nm  = l -  incj
             nm2 = l -2*incj
             np  = l +  incj
#include     "FastS/Compute/SENSOR/wiggle2d_j.for"
#include    "FastC/HPC_LAYER/loopPlan_end.for"
         Endif

      ELSE

#include "FastC/HPC_LAYER/loop_begin.for"                  
            nm  = l -  inck
            nm2 = l -2*inck
            np  = l +  inck
#include    "FastS/Compute/SENSOR/wiggle_k.for"
            nm  = l -  incj
            nm2 = l -2*incj
            np  = l +  incj
#include    "FastS/Compute/SENSOR/wiggle_j.for"
            nm  = l -  inci
            nm2 = l -2*inci
            np  = l +  inci
#include    "FastS/Compute/SENSOR/wiggle_i.for"
c          i = l-lij+ind_loop(1)-1
c        if(k.eq.2.and.j.eq.70.and.ndom.eq.0.and.i.gt.92) then
c         write(*,'(a,7f10.6,i5)')'fluFi',wig(l),test,
c     &  f1,f2,f3,f5, i
c        endif
#include "FastC/HPC_LAYER/loop_end.for"  
                
          if(icorr.eq.1) then !flux manquant en I
            i   = ind_loop(2) + 1
#include    "FastC/HPC_LAYER/loopPlanI_begin.for"
              nm  = l -  inci
              nm2 = l -2*inci
              np  = l +  inci
#include      "FastS/Compute/SENSOR/wiggle_i.for"
#include    "FastC/HPC_LAYER/loopPlan_end.for"
          endif 
          if(jcorr.eq.1) then !complement jmax
           j    = ind_loop(4)+1
#include   "FastC/HPC_LAYER/loopPlanJ_begin.for"                  
              nm  = l -  incj
              nm2 = l -2*incj
              np  = l +  incj
#include      "FastS/Compute/SENSOR/wiggle_j.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"                  
          endif
          if(kcorr.eq.1) then !complement kmax
           k    = ind_loop(6)+1               
#include   "FastC/HPC_LAYER/loopPlanK_begin.for"                  
              nm  = l -  inck
              nm2 = l -2*inck
              np  = l +  inck
#include      "FastS/Compute/SENSOR/wiggle_k.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"                  
          endif

      ENDIF

      end
