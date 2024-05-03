c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine interpolation2(idir, param_int, param_real, ind_loop, 
     &     rop_tmp, rop)
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
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E idir, ind_loop(6), param_int(0:*),ind
      REAL_E param_real(0:*)

      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      REAL_E rop_tmp(param_int(NDIMDX),param_int(NEQ))
    
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq
      REAL_E c1,cv,roe_old,cvinv,ro_old,u_old,v_old,w_old,t_old,roe_tmp
      REAL_E ro_old1,u_old1,v_old1,w_old1,t_old1,roe_old1
 
#include "FastS/formule_param.h"

                          
       neq=param_int(NEQ)
       cv = param_real(CVINF)
       cvinv=1./cv
                   
       !print*, 'cv= ',cv
       !!!!!!! demie-somme de qn et q* !!!!!!!

                  
       do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
         do  i = ind_loop(1), ind_loop(2)   
                              
             l = inddm(i,j,k)
               
             ro_old = rop(l,1)
             u_old  = rop(l,2)
             v_old  = rop(l,3)
             w_old  = rop(l,4)
             t_old  = rop(l,5)

             ro_old1 = rop_tmp(l,1)
             u_old1  = rop_tmp(l,2)
             v_old1  = rop_tmp(l,3)
             w_old1  = rop_tmp(l,4)
             t_old1  = rop_tmp(l,5)

             rop_tmp(l,1)= 0.5*(ro_old + ro_old1)

             rop_tmp(l,2)= 0.5*(ro_old*u_old+ro_old1*u_old1)/
     & rop_tmp(l,1)

             rop_tmp(l,3)= 0.5*(ro_old*v_old+ro_old1*v_old1)/
     & rop_tmp(l,1)

             rop_tmp(l,4)= 0.5*(ro_old*w_old+ro_old1*w_old1)/
     & rop_tmp(l,1)

             roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

             roe_old1 = ro_old1*(cv*t_old1 + 0.5*( u_old1*u_old1
     &                                             +v_old1*v_old1
     &                                             +w_old1*w_old1))



             rop_tmp(l,5)=cvinv*(0.5*(roe_old+roe_old1)/rop_tmp(l,1)-
     &                    0.5*(  rop_tmp(l,2)*rop_tmp(l,2) + 
     &                           rop_tmp(l,3)*rop_tmp(l,3) +
     &                           rop_tmp(l,4)*rop_tmp(l,4) ) )
       
           
                                      
                         
          enddo
         enddo
        enddo    
    

      
                     
    
     
      end
