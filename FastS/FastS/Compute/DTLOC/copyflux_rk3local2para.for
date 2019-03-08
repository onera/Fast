c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copyflux_rk3local2para(param_int,ind_loop,ind_loop_,
     & drodm,stock,ind,taille)
    
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

      INTEGER_E ind_loop(6),ind_loop_(6),param_int(0:*),ind,taille

      REAL_E stock(taille,param_int(NEQ))
      REAL_E drodm(param_int(NDIMDX) ,param_int(NEQ))
    
C Var local

      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3
      REAL_E ratio,coefH,xmut(1),rop(1)
      REAL_E c1,alpha,beta,alphabis

      parameter( alphabis  = -1.2440169358562922 )
      parameter( beta      =  1.8213672050459180 )
      parameter( alpha     = -0.92702963774851155 )

#include "FastS/formule_param.h"
      
      neq=param_int(NEQ)

 
      nistk = (ind_loop_(2)- ind_loop_(1))+1
      nistk2 =(ind_loop_(4)- ind_loop_(3))+1
      nistk3 =(ind_loop_(6)- ind_loop_(5))+1

      if (ind.eq.1) then ! Stockage de f(yn)

       do  ne=1,neq
          do  k = ind_loop(5), ind_loop(6)
           do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)    
                              
              l = inddm(i,j,k)                            
        

               lstk  =  (i+1 - ind_loop_(1))
     &                +(j-ind_loop_(3))*nistk
     &             +(k-ind_loop_(5))*nistk*nistk2


             stock(lstk,ne) =  drodm(l,ne)
                                
           enddo

         enddo
        enddo
       enddo

      elseif (ind.eq.2) then   ! Stockage de alpha*f(yn) + beta*f(y1)
    
       do  ne=1,neq
        do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
          do  i = ind_loop(1), ind_loop(2)                              
               
           l  = inddm(i,j,k)


           lstk  =  (i+1 - ind_loop_(1))
     &              +(j - ind_loop_(3))*nistk
     &            +(k-ind_loop_(5))*nistk*nistk2

           stock(lstk,ne)= (alpha-alphabis)*stock(lstk,ne)
     &                    + beta*drodm(l,ne)
   
           enddo
                          
         enddo
         enddo
         enddo

      end if
      
             
      end
