c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copy_rk3localpara(param_int, ind_loop,ind_loop_,rop, 
     & stock,ind,taille)
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

      INTEGER_E ind_loop(6), ind_loop_(6), param_int(0:*),ind,taille

      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      REAL_E stock(taille,param_int(NEQ))     
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3,cycl
      REAL_E c1

#include "FastS/formule_param.h"

      cycl = param_int(NSSITER)/param_int(LEVEL)                  
      neq=param_int(NEQ)
      
      !print*,"param_int(NEQ)= ",param_int(NEQ) 
      !print*, "ind= ",ind
      !print*, pos
      !print*, taille
      

      nistk  =(ind_loop_(2)- ind_loop_(1))+1
      nistk2 =(ind_loop_(4)- ind_loop_(3))+1
      nistk3 =(ind_loop_(6)- ind_loop_(5))+1

      !print*, ind_loop(1)," ",ind_loop(2)
      !print*, ind_loop(3)," ",ind_loop(4)
      !print*, ind_loop(5)," ",ind_loop(6)

      !print*, "coucou"
       
      if (ind.eq.1) then  ! Stockage 

         
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)

               lstk  =  (i+1 - ind_loop_(1))
     &                +(j-ind_loop_(3))*nistk
     &              +(k-ind_loop_(5))*nistk*nistk2

               stock(lstk,ne) = rop(l,ne)

               !print*, stock(lstk,ne)

              ! print*, 'l= ', l," ",param_int(LEVEL)

        !if (ne==2.and.j==50) then
        !  print*, rop(l,2), ' ',i
        !end if
        !if (j.ge.10.and.ne==1) then
        !   print*, rop(l,ne)
        !end if

      !if(j==ind_loop(4).and.i==ind_loop(2) ) then 
       !print*,"rop stocké= ",stock(lstk,2),k,lstk,ne
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


               lstk  =  (i+1 - ind_loop_(1))
     &                +(j-ind_loop_(3))*nistk
     &              +(k-ind_loop_(5))*nistk*nistk2
                

                  rop(l,ne) =  stock(lstk,ne)

        !if (j==74.and.ne==2.and.i==1) then
        !  print*, "rop= ",rop(l,2)," ",l,""," ",lstk,ind_loop(2)
        !end if

        !if (ne==1.and.j==76.and.i==1.and.param_int(LEVEL)==1) then
        !  print*, rop(l,1)
        !end if


              !  if (j==75.and.ne==1.or.j==74.and.ne==1) then
              !    print*, "ro= ",rop(l,1),"copy"
              ! end if


           enddo

         enddo
         enddo
         enddo

c$$$       if (ind_loop(1)==-1 .or. ind_loop(1)==0) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = ind_loop(3), ind_loop(4)
c$$$                do  i = ind_loop(1), 0
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i+1,j,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$                 
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$
c$$$       else if (ind_loop(2)==param_int(NIJK) .or. 
c$$$     & ind_loop(2)==param_int(NIJK)-1) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = ind_loop(3), ind_loop(4)
c$$$                do  i = param_int(NIJK)-1 , ind_loop(2)
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i-1,j,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$                 
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$          
c$$$
c$$$       else if (ind_loop(3)==-1 .or. ind_loop(3)==0) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = ind_loop(3), 0
c$$$                do  i = ind_loop(1) , ind_loop(2)
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i,j+1,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$                 
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$
c$$$
c$$$       else if (ind_loop(4)==param_int(NIJK+1) .or. 
c$$$     & ind_loop(4)==param_int(NIJK+1)-1) then
c$$$
c$$$         do  ne=1,neq
c$$$           do  k = ind_loop(5), ind_loop(6)
c$$$             do  j = param_int(NIJK+1)-1, ind_loop(4)
c$$$                do  i = ind_loop(1) , ind_loop(2)
c$$$
c$$$                 l  = inddm(i,j,k)
c$$$                 lstk = inddm(i,j-1,k)
c$$$
c$$$                 rop(l,ne) = rop(lstk,ne)
c$$$                 
c$$$                end do
c$$$               end do
c$$$              end do
c$$$            end do
c$$$          
c$$$                   
c$$$        end if
         
    

 
      end if


 

              

      
      end
