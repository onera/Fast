c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine interp_rk3local2(param_int,param_real,coe,ind_loop,
     & stock,stock2,rop,dir,taille,nstep)
   
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

      INTEGER_E ind_loop(6), param_int(0:*),dir,taille,nstep
      INTEGER_E taillefenetre, z(3)
           
      REAL_E param_real(0:*) 
      REAL_E stock(taille,param_int(NEQ))
      REAL_E stock2(taille,param_int(NEQ))
      REAL_E coe(param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk,ind
     &,nistk_
      integer*4 inddm2,i_2,j_2,k_2
      INTEGER_E nistk2,nistk3,l2,lstk2,cycl
      REAL_E c1,roe,cv,cvinv,roe1,roe2
      REAL_E coeff1, coeff2

      
#include "FastS/formule_param.h"  
                        
      ind = param_int(NSSITER)/param_int(LEVEL)
      neq=param_int(NEQ)
       
      cv = param_real(CVINF)
      cvinv=1.0/cv

      !print*, cv

      cycl = param_int(NSSITER)/param_int(LEVEL)


      if (MOD(nstep,cycl)==1) then
         coeff1 = -0.5
         coeff2 = 0.5
      else if (MOD(nstep,cycl)==cycl/2+cycl/4-1) then
         coeff1 = 0.5
         coeff2 = 1.5
      else 
         coeff1 = 0.7886751345948129
         coeff2 = 1.7886751345948129
      end if        

      nistk  = (ind_loop(2) - ind_loop(1)) +1
      nistk2 = (ind_loop(4) - ind_loop(3)) +1
 
      z=0

      if (dir==1 .and. MOD(nstep,cycl).ne.1) then
         z(1)=2
         nistk=4
      else if (dir==2 .and. MOD(nstep,cycl).ne.1) then
         z(2)=2
         nistk2=4
      else if (dir==3 .and. MOD(nstep,cycl).ne.1) then
         z(3)=2
       end if

      if (dir==-1 .and. MOD(nstep,cycl).ne.1) then
         nistk=4
      else if (dir==-2 .and. MOD(nstep,cycl).ne.1) then
         nistk2=4
       end if

      if(dir==1 .and.  MOD(nstep,cycl).ne.1 .and. ind_loop(3).ne.1)then
         
         z(2)=4

      end if


      !print*, z(1),z(2),z(3)

      do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
          do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)


 

               lstk  = (i+1 - ind_loop(1)+z(1))
     &                 +(j- ind_loop(3)+z(2))*nistk
     &                 +(k-ind_loop(5)+z(3))*nistk*nistk2


               lstk2 = (i+1 - ind_loop(1)+z(1))
     &                 +(j- ind_loop(3)+z(2))*nistk
     &                 +(k-ind_loop(5)+z(3))*nistk*nistk2
 

             rop(l,1) = stock2(lstk2,1)*coeff2 -   
     & stock(lstk,1)*coeff1

             rop(l,2) =(stock2(lstk2,1)*stock2(lstk2,2)*coeff2 - 
     & stock(lstk,1)*stock(lstk,2)*coeff1)/rop(l,1)        

             rop(l,3) =(stock2(lstk2,1)*stock2(lstk2,3)*coeff2 - 
     & stock(lstk,1)*stock(lstk,3)*coeff1)/rop(l,1)        

             rop(l,4) =(stock2(lstk2,1)*stock2(lstk2,4)*coeff2 - 
     & stock(lstk,1)*stock(lstk,4)*coeff1)/rop(l,1)    
  

       roe1 = stock(lstk,1)*(cv*stock(lstk,5) + 
     &  0.5*( stock(lstk,2)*stock(lstk,2) + stock(lstk,3)*stock(lstk,3)+
     &   stock(lstk,4)*stock(lstk,4)))

       roe2 = stock2(lstk2,1)*(cv*stock2(lstk2,5) + 
     & 0.5*( stock2(lstk2,2)*stock2(lstk2,2) + 
     & stock2(lstk2,3)*stock2(lstk2,3)+stock2(lstk2,4)*stock2(lstk2,4)))

        rop(l,5) = cvinv*( (roe2*coeff2 - roe1*coeff1)/                              
     &    rop(l,1)  -  0.5*(rop(l,2)*rop(l,2)
     &                + rop(l,3)*rop(l,3)
     &                + rop(l,4)*rop(l,4)))   




      !if(j==75.and.i.ge.147.or.j==76.and.i.ge.147) then 
      ! print*, "lstk= ",lstk,  "lstk2= ",lstk2, "l= ",l
      !end if

              end do                         
            end do
         end do



         end
