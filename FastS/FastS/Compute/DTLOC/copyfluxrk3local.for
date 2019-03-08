c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copyfluxrk3local(idir,param_int,ind_loop,drodm,stock,
     & ind,pos,taille)
    
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

      INTEGER_E idir, ind_loop(6), param_int(0:*),ind,pos,taille

      REAL_E stock(taille/param_int(NEQ),param_int(NEQ))
      REAL_E drodm(param_int(NDIMDX) ,param_int(NEQ))
    
C Var local

      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3
      REAL_E ratio,coefH,xmut(1),rop(1)
      REAL_E c1,alpha,beta,alphabis

      parameter( alphabis =-0.84967936855888582)
      parameter( alpha    =-0.45534180126147938)
      parameter( beta     = 1.2440169358562922)

#include "FastS/formule_param.h"
      
       neq=param_int(NEQ)

 
       nistk = (ind_loop(2)- ind_loop(1))+1
       nistk2 =(ind_loop(4)- ind_loop(3))+1
       nistk3 =(ind_loop(6)- ind_loop(5))+1

      !print*, 'coucouflux',ind       
      
       if (ind.eq.2) then   ! Stockage de alpha*f(yn) + beta*f(y1)
    
      !print*, 'ndom= ', nzone, ' ', 'stockage des flux'      

      do  ne=1,neq
       do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
         do  i = ind_loop(1), ind_loop(2)                              
               
           l  = inddm(i,j,k)


           lstk  =  (i+1 - ind_loop(1))
     &              +(j - ind_loop(3))*nistk
     &            +(k-ind_loop(5))*nistk*nistk2
     &         + pos*nistk*nistk2*nistk3

           stock(lstk,ne) = alpha*stock(lstk,ne) + beta*drodm(l,ne) - 
     & alphabis*stock(lstk,ne)
   
            ! if (lstk.ge.400000) then
            !    print*, lstk, 'nzone= ', nzone
            ! end if

      ! if (j==ind_loop(4).or.j==ind_loop(4)-1.or.
      !&  j==ind_loop(4)-2.or.j==ind_loop(4)-3) then
      !    print*, "stock= ", stock(lstk,1),"  ",j
      !    print*, "drodm= ", drodm(l,1)
      ! end if

           !   if (j==295.and.ne==1) then
           !      print*,stock(lstk,1),"  ",pos,"  ",i
           !   endif
         


           

           enddo
                          
         enddo
         enddo
         enddo


            
      else if (ind.eq.1) then ! Stockage de f(yn)

       !print*, 'ndom= ', nzone, ' ', 'recup des flux'

       do  ne=1,neq
          do  k = ind_loop(5), ind_loop(6)
           do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)    
                              
              l = inddm(i,j,k)                            
        

               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &             +(k-ind_loop(5))*nistk*nistk2
     &         + pos*nistk*nistk2*nistk3


             stock(lstk,ne) =  drodm(l,ne)

         !open(unit=1,file='verifdrodm3.dat')
          !if(ind_loop(1).gt.6) then
          !  write(1,*) drodm(l,ne),"  ",l
          !end if



                                
           enddo

         enddo
        enddo
       enddo

     
       
      end if
      
      !close(1)
             
      end
