c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine interp_rk3local4para(param_int,param_real,coe,ind_loop,
     & ind_loop_,stock,drodmstock,rop,taille,coeff)
   
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

      INTEGER_E ind_loop(6),ind_loop_(6), param_int(0:*),taille
           
      REAL_E param_real(0:*),coeff
      REAL_E stock(taille,param_int(NEQ))
      REAL_E coe(param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      REAL_E drodmstock(taille,param_int(NEQ))
      
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk,ind
      INTEGER_E inddm2,i_2,j_2,k_2
      INTEGER_E nistk2,nistk3,l2,lstk2,posbis
      REAL_E c1,roe,cv,cvinv,ro_old,u_old,v_old,w_old,t_old,roe_old
      REAL_E roinv,u,v,w,dtl

      
#include "FastS/formule_param.h"  
                        
      ind = param_int(NSSITER)/param_int(LEVEL)
      neq=param_int(NEQ)
       
      cv = param_real(CVINF)
      cvinv=1.0/cv

      nistk  = (ind_loop_(2) - ind_loop_(1)) +1
      nistk2 = (ind_loop_(4) - ind_loop_(3)) +1
      nistk3 = (ind_loop_(6) - ind_loop_(5)) +1


      do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
          do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)
 

               lstk  = (i+1 - ind_loop_(1))
     &                 +(j-ind_loop_(3))*nistk
     &                 +(k-ind_loop_(5))*nistk*nistk2


               lstk2 = (i+1 - ind_loop_(1))
     &                 +(j-ind_loop_(3))*nistk
     &                 +(k-ind_loop_(5))*nistk*nistk2

               ro_old = stock(lstk2,1)
               u_old  = stock(lstk2,2)
               v_old  = stock(lstk2,3)
               w_old  = stock(lstk2,4)
               t_old  = stock(lstk2,5)
           

               dtl = coe(l,1)/float(param_int(LEVEL))

              roinv = 1./( stock(lstk2,1)+ coeff*dtl*drodmstock(lstk,1))

         u =(ro_old*stock(lstk2,2) + coeff*dtl*drodmstock(lstk,2))*roinv
         v =(ro_old*stock(lstk2,3) + coeff*dtl*drodmstock(lstk,3))*roinv
         w =(ro_old*stock(lstk2,4) + coeff*dtl*drodmstock(lstk,4))*roinv

         roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                     +v_old*v_old
     &                                     +w_old*w_old))


         rop(l,5)= ( ( roe_old+  coeff*dtl*drodmstock(lstk,5))*roinv
     &              -  0.5*(u*u+v*v+w*w)
     &             )*cvinv   

         rop(l,2) = u
         rop(l,3) = v
         rop(l,4) = w
         rop(l,1) = 1./roinv

              end do                         
            end do
         end do



         end
