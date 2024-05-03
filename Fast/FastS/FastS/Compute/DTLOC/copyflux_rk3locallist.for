c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copyflux_rk3locallist(param_int, drodm, stock, listDnr,
     & nbpts, ind)
    
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

      INTEGER_E param_int(0:*), nbpts, ind
      INTEGER_E listDnr(nbpts)

      REAL_E stock(nbpts,param_int(NEQ))
      REAL_E drodm(param_int(NDIMDX) ,param_int(NEQ))
    
C Var local

      INTEGER_E neq, i, indDnr, ne
      REAL_E alpha,beta,alphabis,alpha_,beta_,alphabis_

      parameter( alphabis   =-0.84967936855888582)
      parameter( alpha      =-0.45534180126147938)
      parameter( beta       = 1.2440169358562922)
      parameter( beta_      =  1.8213672050459180)
      parameter( alpha_     = -0.92702963774851155)
      parameter( alphabis_  = -1.2440169358562922)

#include "FastS/formule_param.h"
      
       neq=param_int(NEQ)

      
      if (ind.eq.2) then   ! Stockage de alpha*f(yn) + beta*f(y1)
    
      !print*, 'ndom= ', nzone, ' ', 'stockage des flux'      

         do  ne=1,neq
            do  i = 1, nbpts 

             indDnr = listDnr(i)                             

             stock(i,ne) = alpha*stock(i,ne) + beta*drodm(indDnr+1,ne) - 
     & alphabis*stock(i,ne)
   
            enddo
         enddo

      else if (ind.eq.3) then   ! Stockage de alpha*f(yn) + beta*f(y1)
    
      !print*, 'ndom= ', nzone, ' ', 'stockage des flux'      

         do  ne=1,neq
            do  i = 1, nbpts 

             indDnr = listDnr(i)                             

             stock(i,ne) = alpha_*stock(i,ne)+ beta_*drodm(indDnr+1,ne)- 
     & alphabis_*stock(i,ne)
   
            enddo
         enddo

            
      else if (ind.eq.1) then ! Stockage de f(yn)

       !print*, 'ndom= ', nzone, ' ', 'recup des flux'

         do  ne=1,neq
            do  i = 1, nbpts 

             indDnr = listDnr(i)            

             stock(i,ne) =  drodm(indDnr+1,ne)

            enddo
         enddo
     
       
      end if
      
      !close(1)
             
      end
