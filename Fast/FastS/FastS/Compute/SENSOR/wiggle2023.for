c***********************************************************************
c     $Date: 2023-07-30 09:20:20 
c     $Revision: 1 $
c     $Author: Ivan Mary $
c***********************************************************************
      subroutine wiggle2023(ndom, param_int, param_real,
     &              ind_grad, ind_sdm, rop, rop_p1, rop_m1, wig, temps)
c***********************************************************************
c
c***********************************************************************
      implicit  none

      REAL_E souszero, cutoff
      parameter(souszero=-1e-12)
      parameter(cutoff=1e-10)

#include "FastS/param_solver.h"


      INTEGER_E ndom,ind_grad(6),ind_sdm(6), param_int(0:*)

      REAL_E rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E rop_m1( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E rop_p1( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E wig( param_int(NDIMDX)*3)
      REAL_E param_real(0:*), temps
    
c Var loc 
      INTEGER_E incmax,i,j,k,l,inci,incj,inck,lij,ltij,lt,lvo,np2,nm3,
     & nm,nm2,np, v1,v2,v3,v4,v5,ind_loop(6),activ,ghost3

      REAL_E qm1,qp1,qm2,qp2,qm3,qp3,qm4,qp4,qm5,qp5,f1,f2,f3,f4,f5,test
      REAL_E f1m,f2m,f3m,f4m,f5m,roref2_inv,vref2_inv,tref2_inv
      REAL_E f1p,f2p,f3p,f4p,f5p
      REAL_E f1n,f2n,f3n,f4n,f5n
      REAL_E f1i,f1im,f1ip
      REAL_E f1i_m,f1im_m,f1ip_m
      REAL_E f1i_p,f1im_p,f1ip_p
      REAL_E sig_m1,sig_m2,sig_p1,sig_p2,h,g

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      inci = 1
      incj = param_int(NIJK)
      inck = param_int(NIJK)*param_int(NIJK+1)

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

      roref2_inv = 1./(param_real(ROINF)*param_real(ROINF))
      vref2_inv  = 1./(param_real(VINF)*param_real(VINF))
      tref2_inv  = 1./(param_real(TINF)*param_real(TINF))

      IF(param_int(ITYPZONE).ne.3)  THEN
       !
       !wig_k
       !
       ind_loop = ind_sdm
       ghost3 = 1
       !gestion bords: 3 point si 2 ghost
       If(param_int(NIJK+4).eq.2) Then
         ghost3 = 0
         if(ind_sdm(5).eq.1) then
           k          = ind_sdm(5)
           ind_loop(5)= 2
#include   "FastC/HPC_LAYER/loopPlanK_begin.for"
              nm  = l -  inck
              nm2 = l -2*inck
              np  = l +  inck
#include      "FastS/Compute/SENSOR/wiggle_k.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"
         endif
         if(ind_sdm(6).eq.param_int(IJKV+2)) then
           k   =  ind_sdm(6)+1
#include   "FastC/HPC_LAYER/loopPlanK_begin.for"
              nm  = l -  inck
              nm2 = l -2*inck
              np  = l +  inck
#include      "FastS/Compute/SENSOR/wiggle_k.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"
         endif
       Endif
       !gestion coeur domaine k
       if(ind_sdm(6).eq.param_int(IJKV+2).and.ghost3.eq.1)
     &      ind_loop(6)=ind_sdm(6)+1
#include "FastC/HPC_LAYER/loop_begin.for"                  
            nm  = l -  inck
            nm3 = l -3*inck
            nm2 = l -2*inck
            np  = l +  inck
            np2 = l +2*inck
#include    "FastS/Compute/SENSOR/wiggle5_k.for"
#include "FastC/HPC_LAYER/loop_end.for"                  
      Endif !zone 3d

       !
       !wig_j
       !
       ind_loop = ind_sdm
       ghost3 = 1
       !gestion bords j: 3 point si 2 ghost
       If(param_int(NIJK+3).eq.2) Then
         ghost3 = 0
         if(ind_sdm(3).eq.1) then
           j          =  ind_sdm(3)
           ind_loop(3)= 2
#include   "FastC/HPC_LAYER/loopPlanJ_begin.for" 
              nm  = l -  incj
              nm2 = l -2*incj
              np  = l +  incj
#include      "FastS/Compute/SENSOR/wiggle_j.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"
         endif
         if(ind_sdm(4).eq.param_int(IJKV+1)) then
           j =  ind_sdm(4)+1
#include   "FastC/HPC_LAYER/loopPlanJ_begin.for" 
              nm  = l -  incj
              nm2 = l -2*incj
              np  = l +  incj
#include    "FastS/Compute/SENSOR/wiggle_j.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"
         endif
       Endif
       !gestion coeur domaine j
       if(ind_sdm(4).eq.param_int(IJKV+1).and.ghost3.eq.1)
     &  ind_loop(4)=ind_sdm(4)+1
#include "FastC/HPC_LAYER/loop_begin.for"                  
            nm  = l -  incj
            nm3 = l -3*incj
            nm2 = l -2*incj
            np  = l +  incj
            np2 = l +2*incj
#include    "FastS/Compute/SENSOR/wiggle5_j.for"
#include "FastC/HPC_LAYER/loop_end.for"                  
       !
       !wig_i
       !
       ind_loop = ind_sdm
       ghost3 = 1
       !gestion bords i: 3 point si 2 ghost
       If(param_int(NIJK+3).eq.2) Then
         ghost3 = 0
         if(ind_sdm(1).eq.1) then
           i          =  ind_sdm(1)
           ind_loop(1)= 2
#include   "FastC/HPC_LAYER/loopPlanI_begin.for" 
              nm  = l -  inci
              nm2 = l -2*inci
              np  = l +  inci
#include      "FastS/Compute/SENSOR/wiggle_i.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"
         endif
         if(ind_sdm(2).eq.param_int(IJKV)) then
           i =  ind_sdm(2)+1
#include   "FastC/HPC_LAYER/loopPlanI_begin.for" 
              nm  = l -  inci
              nm2 = l -2*inci
              np  = l +  inci
#include      "FastS/Compute/SENSOR/wiggle_i.for"
#include   "FastC/HPC_LAYER/loopPlan_end.for"
         endif
       Endif
       !gestion coeur domaine i
       if(ind_sdm(2).eq.param_int(IJKV).and.ghost3.eq.1)
     &     ind_loop(2)=ind_sdm(2)+1

#include "FastC/HPC_LAYER/loop_begin.for"                  
          nm  = l -  inci
          nm3 = l -3*inci
          nm2 = l -2*inci
          np  = l +  inci
          np2 = l +2*inci
#include  "FastS/Compute/SENSOR/wiggle5_i.for"
#include "FastC/HPC_LAYER/loop_end.for"                  

      end
