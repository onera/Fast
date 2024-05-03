c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine switchvectors(rop,roptmp,param_int)
    
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    extrapolation ordere zero cell fictive
c
c     VAL
c_V    Optimisation NEC
c
c     COM
c***********************************************************************

!#include "parallelF.h"     
      implicit none


#include "FastS/param_solver.h"

      INTEGER_E param_int(0:*), indloop(6),i,j,k,l,ne,neq,cycl

      REAL_E rop(param_int(NDIMDX) ,param_int(NEQ))
      REAL_E roptmp(param_int(NDIMDX) ,param_int(NEQ))
       !REAL_E ropm1(param_int(NDIMDX) ,param_int(NEQ))
      REAL_E stk
    

      
      !print*, rop

#include "FastS/formule_param.h"

      indloop(1)=1
      indloop(2)=param_int(IJKV)
      indloop(3)=1
      indloop(4)=param_int(IJKV+1)
      indloop(5)=1
      indloop(6)=param_int(IJKV+2)
                        
      neq=param_int(NEQ)

      cycl=param_int(NSSITER)/param_int(LEVEL)


      !print*,indloop(2)
      !print*,indloop(4)
      !print*,indloop(6)

      !print*, 'coucou', param_int(LEVEL)

      do  ne=1,neq
        do  k = indloop(5), indloop(6)
          do  j = indloop(3), indloop(4)
            do  i = indloop(1), indloop(2)   
                              
               l  = inddm(i,j,k)
               
               !print*, l

               stk = rop(l,ne)
               rop(l,ne) = roptmp(l,ne)
               roptmp(l,ne) = stk


       ! if (j==75.and.i==indloop(1).and.cycl==16.and.ne==1) then
       !    print*, rop(l,1)
       ! end if
               
             end do
            end do
          end do
        end do
         
        !print*, "coucou",param_int(LEVEL) 





      end





