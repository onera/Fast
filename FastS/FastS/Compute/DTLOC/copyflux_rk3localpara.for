c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copyflux_rk3localpara(param_int,ind_loop,ind_loop_,
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

      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,neq,n,nistk,lstk,lvec
      INTEGER_E nistk2,nistk3
      REAL_E ratio,coefH,xmut(1),rop(1)
      REAL_E c1,alpha,beta,alphabis, tmp(6400,param_int(NEQ))

      parameter( alphabis =-0.84967936855888582)
      parameter( alpha    =-0.45534180126147938)
      parameter( beta     = 1.2440169358562922)

#include "FastS/formule_param.h"
      
       neq=param_int(NEQ)

 
       nistk = (ind_loop_(2)- ind_loop_(1))+1
       nistk2 =(ind_loop_(4)- ind_loop_(3))+1
       nistk3 =(ind_loop_(6)- ind_loop_(5))+1

      
       if (ind.eq.2) then   ! Stockage de alpha*f(yn) + beta*f(y1)
    
       if(neq.eq.5) then

         do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)
               
             lstk  =  (1+i - ind_loop_(1))
     &               +(  j - ind_loop_(3))*nistk
     &               +(  k - ind_loop_(5))*nistk*nistk2

             lvec= 1+i - ind_loop_(1)

             tmp(lvec,1) = (alpha-alphabis)*stock(lstk,1)
             tmp(lvec,2) = (alpha-alphabis)*stock(lstk,2)
             tmp(lvec,3) = (alpha-alphabis)*stock(lstk,3)
             tmp(lvec,4) = (alpha-alphabis)*stock(lstk,4)
             tmp(lvec,5) = (alpha-alphabis)*stock(lstk,5)
            enddo
            do  i = ind_loop(1), ind_loop(2)
               
             l     = inddm(i,j,k)
             lvec  = 1+i - ind_loop_(1)

             lstk  =  (1+i - ind_loop_(1))
     &               +(  j - ind_loop_(3))*nistk
     &               +(  k - ind_loop_(5))*nistk*nistk2

             stock(lstk,1) = tmp(lvec,1) + beta*drodm(l,1)
             stock(lstk,2) = tmp(lvec,2) + beta*drodm(l,2)
             stock(lstk,3) = tmp(lvec,3) + beta*drodm(l,3)
             stock(lstk,4) = tmp(lvec,4) + beta*drodm(l,4)
             stock(lstk,5) = tmp(lvec,5) + beta*drodm(l,5)
            enddo
         enddo
        enddo

       else

         do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)
               
             lstk  =  (1+i - ind_loop_(1))
     &               +(  j - ind_loop_(3))*nistk
     &               +(  k - ind_loop_(5))*nistk*nistk2

             lvec= 1+i - ind_loop_(1)

             tmp(lvec,1) = (alpha-alphabis)*stock(lstk,1)
             tmp(lvec,2) = (alpha-alphabis)*stock(lstk,2)
             tmp(lvec,3) = (alpha-alphabis)*stock(lstk,3)
             tmp(lvec,4) = (alpha-alphabis)*stock(lstk,4)
             tmp(lvec,5) = (alpha-alphabis)*stock(lstk,5)
             tmp(lvec,6) = (alpha-alphabis)*stock(lstk,6)
            enddo
            do  i = ind_loop(1), ind_loop(2)
               
             l     = inddm(i,j,k)
             lvec  = 1+i - ind_loop_(1)
 
             lstk  =  (1+i - ind_loop_(1))
     &               +(  j - ind_loop_(3))*nistk
     &               +(  k - ind_loop_(5))*nistk*nistk2

             stock(lstk,1) = tmp(lvec,1) + beta*drodm(l,1)
             stock(lstk,2) = tmp(lvec,2) + beta*drodm(l,2)
             stock(lstk,3) = tmp(lvec,3) + beta*drodm(l,3)
             stock(lstk,4) = tmp(lvec,4) + beta*drodm(l,4)
             stock(lstk,5) = tmp(lvec,5) + beta*drodm(l,5)
             stock(lstk,6) = tmp(lvec,6) + beta*drodm(l,6)
            enddo
          enddo
         enddo

       endif


            
      else if (ind.eq.1) then ! Stockage de f(yn)

       !print*, 'ndom= ', nzone, ' ', 'recup des flux'

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

     
       
      end if
      
      !close(1)
             
      end
