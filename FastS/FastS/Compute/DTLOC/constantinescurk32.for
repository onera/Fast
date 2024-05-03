c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine constantinescurk32(idir,param_int,param_real,ind_loop, 
     &     rop,iptstk,nzone)
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

      INTEGER_E idir, ind_loop(6), param_int(0:*),nzone
      REAL_E param_real(0:*)

      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      REAL_E iptstk(4000000,param_int(NEQ))
      REAL_E stock(200000,param_int(NEQ))
      
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3,lstk2
      REAL_E c1,cv,rhoe,rhoe_tmp,cvinv,ro_old,u_old,v_old,w_old,t_old,
     & roe_old,roe_stk

#include "FastS/formule_param.h"

                          
       neq=param_int(NEQ)
       cv = param_real(CVINF)
       cvinv=1./cv
      
            
       nistk  = (ind_loop(2)- ind_loop(1))+1
       nistk2 = (ind_loop(4)- ind_loop(3))+1
       nistk3 = (ind_loop(6)- ind_loop(5))+1
        

       !!!!!!! stockage de rop !!!!!!!

                  
        ! do  ne=1,neq
           !  do  k = ind_loop(5), ind_loop(6)
           !    do  j = ind_loop(3), ind_loop(4)
           !     do  i = ind_loop(1), ind_loop(2)   
                              
           !    l  = inddm(i,j,k)

          !     lstk  =  (i+1 - ind_loop(1))
     &    !            +(j-ind_loop(3))*nistk
     &    !          +(k-ind_loop(5))*nistk*nistk2

         !      stock(lstk,ne) = rop(l,ne)
                                      
        !   enddo
                          
        ! enddo
        ! enddo
        !enddo    
    
       

        !! Modif des 4 premieres ou dernieres cellules !!

        
        do  k = ind_loop(5), ind_loop(6)
           do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)   
                 
                 l  = inddm(i,j,k)

                 lstk2  =  (i+1 - ind_loop(1))
     &                  +(j-ind_loop(3))*nistk
     &                  +(k-ind_loop(5))*nistk*nistk2
     &                  + nzone*nistk*nistk2*nistk3


               ro_old = rop(l,1)
               u_old  = rop(l,2)
               v_old  = rop(l,3)
               w_old  = rop(l,4)
               t_old  = rop(l,5)

               rop(l,1)=0.5*(rop(l,1)+iptstk(lstk2,1))

               rop(l,2)=0.5*(ro_old*rop(l,2)+
     &               iptstk(lstk2,1)*iptstk(lstk2,2))/rop(l,1)

               rop(l,3)=0.5*(ro_old*rop(l,3)+
     &               iptstk(lstk2,1)*iptstk(lstk2,3))/rop(l,1)

               rop(l,4)=0.5*(ro_old*rop(l,4)+
     &               iptstk(lstk2,1)*iptstk(lstk2,4))/rop(l,1)

             roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

                roe_stk = iptstk(lstk2,1)*(cv*iptstk(lstk2,5)
     &                  + 0.5*(iptstk(lstk2,2)*iptstk(lstk2,2)
     &                  + iptstk(lstk2,3)*iptstk(lstk2,3)
     &                  + iptstk(lstk2,4)*iptstk(lstk2,4) ) )


               rop(l,5)=cvinv*(0.5*(roe_old+roe_stk)/rop(l,1)-
     &                    0.5*(  rop(l,2)*rop(l,2) + 
     &                           rop(l,3)*rop(l,3) +
     &                           rop(l,4)*rop(l,4) ) )


              ! if (j==50) then
              !    print*, iptstk(lstk2,5), rop(l,5), t_old
              ! end if

 
              enddo
              
           enddo
        enddo                    
    
     
      end
