c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copy(idir, param_int, ind_loop, rop, stock, ind,nzone)
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

      INTEGER_E idir, ind_loop(6), param_int(0:*),ind,nzone

      REAL_E rop( param_int(NDIMDX),param_int(NEQ))
      REAL_E stock(4000000,param_int(NEQ))
      
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3
      REAL_E c1

#include "FastS/formule_param.h"

                          
       neq=param_int(NEQ)


       nistk = (ind_loop(2)- ind_loop(1))+1
       nistk2 =(ind_loop(4)- ind_loop(3))+1
       nistk3 =(ind_loop(6)- ind_loop(5))+1
       
       if (ind.eq.1) then  ! Stockage 
         
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)

               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3

               stock(lstk,ne) = rop(l,ne)
               !print*, stock(lstk,1)-rop(l,1)

               !if (j==50.and.ne==2) then
               !   print*, stock(lstk,ne), rop(l,ne)
               !end if
                         
           enddo
                          
         enddo
         enddo
         enddo
         

      elseif (ind.eq.2) then  ! Recuperation des valeurs
      
         do  ne=1,neq
           do  k = ind_loop(5), ind_loop(6)
             do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)               
               
                  l  = inddm(i,j,k)


               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3

                

                  rop(l,ne) =  stock(lstk,ne)
                                   
           enddo

         enddo
         enddo
         enddo
         
    

 
      end if

              

      
      end
