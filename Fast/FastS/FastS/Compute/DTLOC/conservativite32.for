c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine conservativite32(idir,param_int,param_real,param_int2, 
     & ind_loop,imax,rop,drodmzone,drodm,stock,coe,nstep,nzone)
    
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

      INTEGER_E idir, ind_loop(6), param_int(0:*),nstep,nzone,imax
      INTEGER_E param_int2(0:*)
      
      REAL_E param_real(0:*)
      REAL_E rop( param_int(NDIMDX),param_int(NEQ))
      REAL_E drodmzone(param_int(NDIMDX),param_int(NEQ))
      REAL_E drodm(param_int2(NDIMDX),param_int2(NEQ))
      REAL_E stock(4000000,param_int(NEQ))
      REAL_E coe(param_int(NDIMDX),param_int(NEQ_COE))
     
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk,ind
      INTEGER_E inddm2,i_2,j_2,k_2
      INTEGER_E nistk2,nistk3,l2,lstk1
      REAL_E c1,roe,cv,cvinv,ro_old,u_old,v_old,w_old,t_old,roe_old

        inddm2(i_2,j_2,k_2) = 1  
     &  + (i_2+param_int2(NIJK+3)-1)
     &  + (j_2+param_int2(NIJK+3)-1)*param_int2(NIJK)
     &  + (k_2+param_int2(NIJK+4)-1)*param_int2(NIJK)*param_int2(NIJK+1)
      
#include "FastS/formule_param.h"                          
       ind = param_int(NSSITER)/param_int(LEVEL)
       neq=param_int(NEQ)
       
       cv = param_real(CVINF)
       cvinv=1.0/cv
            

       nistk  = (ind_loop(2) - ind_loop(1)) +1
       nistk2 = (ind_loop(4) - ind_loop(3)) +1
       nistk3 = (ind_loop(6) - ind_loop(5)) +1
       
       if (mod(nstep,ind)==ind/2) then  !!! Recuperation de Fi(**)-Fi(n)
         
          if (ind_loop(1).ne.1) then
           do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)   
                              
                   !l=inddm(i,j,k)
                   l2=inddm2(imax,j,k)

               lstk  =  (i+1 - ind_loop(1))
     &              +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &              + nzone*nistk*nistk2*nistk3
     
        
               stock(lstk,ne)=drodm(l2-1,ne)

            
            
    

           enddo
                          
         enddo
         enddo
         enddo
        
          else if(ind_loop(1).eq.1) then   !!! Recuperation de Fi(**)-Fi(n)
           
            do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)   
                              
               !l  = inddm(i,j,k)
               l2=inddm2(imax,j,k)
               

               lstk  =  (i+1 - ind_loop(1))
     &                +(j- ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3
               
                          
               stock(lstk,ne)=drodm(l2+1,ne)

 
                          
                 
           enddo
                          
         enddo
         enddo
         enddo
         
        end if

      end if 
 
       
        if (mod(nstep,ind)==0) then
        
    
      !!! Modif de la derniere maille de la zone !!!


         if (ind_loop(1).ne.1) then
       
           do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)
               l2 = inddm2(imax,j,k)

               lstk  = (i+1 - ind_loop(1))
     &                 +(j- ind_loop(3))*nistk
     &                 +(k-ind_loop(5))*nistk*nistk2
     &                 + nzone*nistk*nistk2*nistk3


              ro_old = rop(l,1)
              u_old  = rop(l,2)
              v_old  = rop(l,3)
              w_old  = rop(l,4)
              t_old  = rop(l,5)
           

               rop(l,1) = rop(l,1) 
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2-1,1)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,1)*0.25


               rop(l,2) =(ro_old*rop(l,2)
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2-1,2)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,2)*0.25)/rop(l,1)

       

               rop(l,3) =(ro_old*rop(l,3)
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2-1,3)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,3)*0.25)/rop(l,1)
       

               rop(l,4) =(ro_old*rop(l,4)
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2-1,4)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,4)*0.25)/rop(l,1)



             roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))


              rop(l,5) = cvinv*( ( roe_old -
     &(coe(l,1)/float(param_int(LEVEL)))*drodm(l2-1,5)*0.25 +
     &(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,5)*0.25 )/rop(l,1)                              
     &  -  0.5*(rop(l,2)*rop(l,2)
     &        + rop(l,3)*rop(l,3)
     &        + rop(l,4)*rop(l,4)))   



         
              enddo                         
            enddo
         enddo

       end if
        
      
       !!! Modif de la premiere maille de la zone !!!


      if (ind_loop(1).eq.1) then

           do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)
               l2 = inddm2(imax,j,k)

               lstk  = (i+1 - ind_loop(1))
     &                 +(j- ind_loop(3))*nistk
     &                 +(k-ind_loop(5))*nistk*nistk2
     &                 + nzone*nistk*nistk2*nistk3


              ro_old = rop(l,1)
              u_old  = rop(l,2)
              v_old  = rop(l,3)
              w_old  = rop(l,4)
              t_old  = rop(l,5)
           

               rop(l,1) = rop(l,1) 
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2+1,1)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,1)*0.25

 
               rop(l,2) =(ro_old*rop(l,2)
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2+1,2)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,2)*0.25)/rop(l,1)        

               rop(l,3) =(ro_old*rop(l,3)
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2+1,3)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,3)*0.25)/rop(l,1)        

               rop(l,4) =(ro_old*rop(l,4)
     & -(coe(l,1)/float(param_int(LEVEL)))*drodm(l2+1,4)*0.25
     & +(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,4)*0.25)/rop(l,1)  

             roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))


              rop(l,5) = cvinv*( ( roe_old -
     &(coe(l,1)/float(param_int(LEVEL)))*drodm(l2+1,5)*0.25 +
     &(coe(l,1)/float(param_int(LEVEL)))*stock(lstk,5)*0.25 )/rop(l,1)                              
     &  -  0.5*(rop(l,2)*rop(l,2)
     &        + rop(l,3)*rop(l,3)
     &        + rop(l,4)*rop(l,4)))   

         
              enddo                         
            enddo
         enddo

       end if
 
        
      end if    
 
                 

      
      end
