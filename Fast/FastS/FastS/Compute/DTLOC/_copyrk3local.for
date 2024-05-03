c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copy_rk3local(param_int, ind_loop, rop, stock, 
     & ind,taille)
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

      INTEGER_E ind_loop(6), param_int(0:*),ind,taille

      REAL_E rop( param_int(NDIMDX),param_int(NEQ))
      REAL_E stock(taille/param_int(NEQ),param_int(NEQ))     
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3,cycl
      REAL_E c1

#include "FastS/formule_param.h"

      cycl = param_int(NSSITER)/param_int(LEVEL)                  
      neq=param_int(NEQ)
      
      !print*, neq
      !print*, "ind= ",ind
      !print*, pos
      !print*, taille
      

      nistk  =(ind_loop(2)- ind_loop(1))+1
      nistk2 =(ind_loop(4)- ind_loop(3))+1
      nistk3 =(ind_loop(6)- ind_loop(5))+1

      !print*, ind_loop(1)," ",ind_loop(2)
      !print*, ind_loop(3)," ",ind_loop(4)
      !print*, ind_loop(5)," ",ind_loop(6)
       
      if (ind.eq.1) then  ! Stockage 

         
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)

               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2

               stock(lstk,ne) = rop(l,ne)

        !if (j.ge.10.and.ne==1) then
        !   print*, rop(l,ne)
        !end if
                          
           enddo
                          
         enddo
         enddo
         enddo
         

      elseif (ind.eq.2) then  ! Recuperation des valeurs

         !print*, "coucou recup valeurs"
      
         do  ne=1,neq
           do  k = ind_loop(5), ind_loop(6)
             do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)               
               
                  l  = inddm(i,j,k)


               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
                

                  rop(l,ne) =  stock(lstk,ne)

        !if (j==45.and.ne==1) then
        !  print*, "rop= ",rop(l,1)," ",l," ",lstk," ",ind," ",pos
        !end if


               ! if (j==45.and.ne==1) then
               !   print*, rop(l,1)," ",l,"  ",lstk
               !end if


           enddo

         enddo
         enddo
         enddo
         
    

 
      end if

              

      
      end
