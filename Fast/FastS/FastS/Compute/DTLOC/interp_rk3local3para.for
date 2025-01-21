c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine interp_rk3local3para( param_int,param_real,coe,
     &                                 ind_loop, ind_loop_,
     &                                 stock, stock2, rop,
     &                                 dir,taille, nstep)
   
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

      INTEGER_E ind_loop(6),ind_loop_(6),param_int(0:*),dir,taille,nstep
      INTEGER_E taillefenetre
           
      REAL_E param_real(0:*) 
      REAL_E stock(taille,param_int(NEQ))
      REAL_E stock2(taille,param_int(NEQ))
      REAL_E coe(param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
     &,nistk_,  z(3)
      INTEGER_E inddm2,i_2,j_2,k_2
      INTEGER_E nistk2,nistk3,l2,lstk2,cycl
      REAL_E c1,roe,cv,cvinv,roe1,roe2
      REAL_E coeff1, coeff2,roinv,u,v,w

      
#include "FastS/formule_param.h"  
                        
      neq  = param_int(NEQ)
      cycl = param_int(NSSITER)/param_int(LEVEL)
       
      cv = param_real(CVINF)
      cvinv=1.0/cv



      if (MOD(nstep,cycl)==1) then
         coeff1 = -0.5
         coeff2 = 0.5
      else if (MOD(nstep,cycl)==cycl/4) then
         coeff1 = -0.21132486540518713
         coeff2 = 0.78867513459481287
      else if (MOD(nstep,cycl)==cycl/2+1) then
         coeff1 = 0.5
         coeff2 = 1.5
      end if        

      nistk  = (ind_loop_(2) - ind_loop_(1)) +1
      nistk2 = (ind_loop_(4) - ind_loop_(3)) +1
 
      z=0

      do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
          do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)

               lstk  = (1+i - ind_loop_(1)+z(1))
     &                +(  j - ind_loop_(3)+z(2))*nistk
     &                +(  k  -ind_loop_(5)+z(3))*nistk*nistk2


               lstk2 = (1+i- ind_loop_(1)+z(1))
     &                +(  j- ind_loop_(3)+z(2))*nistk
     &                +(  k- ind_loop_(5)+z(3))*nistk*nistk2
 

             roinv = 1./(stock2(lstk2,1)*coeff2 -  stock(lstk,1)*coeff1)

             u = ( stock2(lstk2,1)*stock2(lstk2,2)*coeff2 
     &            -stock(  lstk,1)*stock( lstk ,2)*coeff1 )*roinv

             v = ( stock2(lstk2,1)*stock2(lstk2,3)*coeff2 
     &            -stock(  lstk,1)*stock(  lstk,3)*coeff1)*roinv

             w = ( stock2(lstk2,1)*stock2(lstk2,4)*coeff2 
     &            -stock(  lstk,1)*stock(  lstk,4)*coeff1)*roinv


             roe1 = stock(lstk,1)*(     cv*stock(lstk,5)  
     &                             +0.5*( stock(lstk,2)*stock(lstk,2) 
     &                                   +stock(lstk,3)*stock(lstk,3)
     &                                   +stock(lstk,4)*stock(lstk,4) )
     &                            )

            roe2=stock2(lstk2,1)*(     cv*stock2(lstk2,5) 
     &                           +0.5*( stock2(lstk2,2)*stock2(lstk2,2)
     &                                 +stock2(lstk2,3)*stock2(lstk2,3)
     &                                 +stock2(lstk2,4)*stock2(lstk2,4))
     &                           )

            rop(l,5) = cvinv*( (roe2*coeff2 - roe1*coeff1)*roinv
     &                        -  0.5*( u*u+v*v+w*w)
     &                       )   

            rop(l,1) = 1./roinv
            rop(l,2) = u
            rop(l,3) = v
            rop(l,4) = w

              end do                         
            end do
         end do



         end
