c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine interp_rk3local3list(param_int,param_real,listRcv,
     & stock,stock2,rop,nbpts,nstep)
   
c***********************************************************************
c_U   USER : PECHIER ********SUBROUTINE FONCTIONNANT SEULE**************
c
c     ACT
c_A    extrapolation ordere zero cell fictive
c
c     VAL
c_V    Optimisation NEC
c
c     COM
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E param_int(0:*),nbpts,nstep
      INTEGER_E listRcv(nbpts)
           
      REAL_E param_real(0:*) 
      REAL_E stock(nbpts,param_int(NEQ))
      REAL_E stock2(nbpts,param_int(NEQ))
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      
 
C Var local

      INTEGER_E neq,i,indDnr,ne,cycl,indRcv
      REAL_E coeff1, coeff2, cv, cvinv,roe1,roe2

      
#include "FastS/formule_param.h"  
                        
      neq=param_int(NEQ)
       
      cv = param_real(CVINF)
      cvinv=1.0/cv

      !print*, cv

      cycl = param_int(NSSITER)/param_int(LEVEL)


      if (MOD(nstep,cycl)==1) then
         coeff1 = -0.5
         coeff2 = 0.5
      else if (MOD(nstep,cycl)==cycl/4) then
         coeff1 = -0.21132486540518713
         coeff2 = 0.78867513459481287
      else if (MOD(nstep,cycl)==cycl/2+cycl/4-1) then
         coeff1 = 0.5
         coeff2 = 1.5
      end if        

 

      do i=1,nbpts

         indRcv = listRcv(i)

 
            rop(indRcv+1,1) = stock2(i,1)*coeff2 -   
     & stock(i,1)*coeff1

            rop(indRcv+1,2) =(stock2(i,1)*stock2(i,2)*coeff2 - 
     & stock(i,1)*stock(i,2)*coeff1)/rop(indRcv+1,1)        

            rop(indRcv+1,3) =(stock2(i,1)*stock2(i,3)*coeff2 - 
     & stock(i,1)*stock(i,3)*coeff1)/rop(indRcv+1,1)        

            rop(indRcv+1,4) =(stock2(i,1)*stock2(i,4)*coeff2 - 
     & stock(i,1)*stock(i,4)*coeff1)/rop(indRcv+1,1)    
  

       roe1 = stock(i,1)*(cv*stock(i,5) + 
     &  0.5*( stock(i,2)*stock(i,2) + stock(i,3)*stock(i,3)+
     &   stock(i,4)*stock(i,4)))

       roe2 = stock2(i,1)*(cv*stock2(i,5) + 
     & 0.5*( stock2(i,2)*stock2(i,2) + 
     & stock2(i,3)*stock2(i,3)+stock2(i,4)*stock2(i,4)))

        rop(indRcv+1,5) = cvinv*( (roe2*coeff2 - roe1*coeff1)/                              
     &    rop(indRcv+1,1) - 0.5*(rop(indRcv+1,2)*rop(indRcv+1,2)
     &                + rop(indRcv+1,3)*rop(indRcv+1,3)
     &                + rop(indRcv+1,4)*rop(indRcv+1,4)))   




      end do



       end
