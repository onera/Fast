c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine constantinescurk3(idir,param_int,param_real,ind_loop, 
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
      REAL_E iptstk(400000,param_int(NEQ))
      REAL_E stock(200000,param_int(NEQ))
      
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3,lstk2
      REAL_E c1,cv,rhoe,rhoe_tmp,cvinv,a

#include "FastS/formule_param.h"

                          
       neq=param_int(NEQ)
       cv = param_real(CVINF)
       cvinv=1./cv
       a=1./2.
            
       nistk  = (ind_loop(2)- ind_loop(1))+1
       nistk2 = (ind_loop(4)- ind_loop(3))+1
       nistk3 = (ind_loop(6)- ind_loop(5))+1
        

       !!!!!!! stockage de rop !!!!!!!

                  
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                do  i = ind_loop(1), ind_loop(2)   
                              
               l  = inddm(i,j,k)

               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2

               stock(lstk,ne) = rop(l,ne)
                                      
           enddo
                          
         enddo
         enddo
        enddo    
    
       

        !!!!!!! Densite !!!!!!

        
        do  k = ind_loop(5), ind_loop(6)
           do  j = ind_loop(3), ind_loop(4)
              do  i = ind_loop(1), ind_loop(2)   
                 
                 l  = inddm(i,j,k)

               lstk2  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3



                 rop(l,1)=a*(rop(l,1)+iptstk(lstk2,1))

 
              enddo
              
           enddo
        enddo

        !!!!!!!! u,v,w !!!!!!

        do ne=2,4
           do  k = ind_loop(5), ind_loop(6)
              do  j = ind_loop(3), ind_loop(4)
                 do  i = ind_loop(1), ind_loop(2)   
                    
                    l  = inddm(i,j,k)

                    lstk  =  (i+1 - ind_loop(1))
     &                   +(j-ind_loop(3))*nistk
     &              + (k-ind_loop(5))*nistk*nistk2

               lstk2  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3



                    rop(l,ne)=a*(stock(lstk,1)*stock(lstk,ne)+
     &               iptstk(lstk2,1)*iptstk(lstk2,ne) )/rop(l,1)


                    
                 enddo
                 
              enddo
           enddo
          enddo


         !!!!!!! T !!!!!!

           
           do  k = ind_loop(5), ind_loop(6)
              do  j = ind_loop(3), ind_loop(4)
                 do  i = ind_loop(1), ind_loop(2)   
                    
                    l  = inddm(i,j,k)

                    lstk  =  (i+1 - ind_loop(1))
     &                   +(j-ind_loop(3))*nistk
     &              + (k-ind_loop(5))*nistk*nistk2


               lstk2  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &              +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3


          ! if (j==55.and.i==ind_loop(2)) then
          !    print*, ne," y6= ", rop(l,5)," ","y3= ",iptstk(lstk2,5) 
          ! end if     

                    rhoe_tmp= stock(lstk,1)*( cv*stock(lstk,5)
     &                   + a*( stock(lstk,2)*stock(lstk,2)
     &                   + stock(lstk,3)*stock(lstk,3)
     &                   + stock(lstk,4)*stock(lstk,4) ) )

                    rhoe = iptstk(lstk2,1)*(cv*iptstk(lstk2,5)
     &                   + a*(iptstk(lstk2,2)*iptstk(lstk2,2)
     &                   + iptstk(lstk2,3)*iptstk(lstk2,3)
     &                   + iptstk(lstk2,4)*iptstk(lstk2,4) ) )

                rop(l,5)=cvinv*(a*(rhoe_tmp+rhoe)/rop(l,1)-
     &                    a*(  rop(l,2)*rop(l,2) + 
     &                         rop(l,3)*rop(l,3) +
     &                         rop(l,4)*rop(l,4) ) )


          ! if (j==55.and.i==ind_loop(2)) then
          !    print*, ne," yn+1= ", rop(l,5) 
          ! end if     
                    
                 enddo
                 
              enddo
           enddo


          ! do  k = ind_loop(5), ind_loop(6)
          !    j=50
          !      do  i = ind_loop(1)-2, ind_loop(2)   
                              
           !    l  = inddm(i,j,k)

            !   print*, rop_tmp(l,5),i 
                                      
          ! enddo
                          
        ! enddo
         
    

                    
    
     
      end
