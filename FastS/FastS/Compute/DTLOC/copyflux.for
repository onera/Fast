c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copyflux(idir,param_int,ind_loop,drodm,stock,ind,nzone)
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

      INTEGER_E idir, ind_loop(6), param_int(0:*),ind,nzone

      REAL_E stock(4000000 ,param_int(NEQ))
      REAL_E drodm(param_int(NDIMDX) ,param_int(NEQ))
    
C Var local

      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,neq,n,nistk,lstk
      INTEGER_E nistk2,nistk3
      REAL_E ratio,coefH,xmut(1),rop(1)
      REAL_E c1

#include "FastS/formule_param.h"
      
            neq=param_int(NEQ)

       nistk = (ind_loop(2)- ind_loop(1))+1
       nistk2 =(ind_loop(4)- ind_loop(3))+1
       nistk3 =(ind_loop(6)- ind_loop(5))+1

      !print*, 'coucouflux',ind       
      

      if (ind.eq.1) then   ! Stockage
    
            

      do  ne=1,neq
       do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
         do  i = ind_loop(1), ind_loop(2)                              
               
              l  = inddm(i,j,k)


               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &             +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3

              stock(lstk,ne) = drodm(l,ne)


           enddo
                          
         enddo
         enddo
         enddo


            
      else if (ind.eq.2) then ! recuperation des valeurs

       do  ne=1,neq
          do  k = ind_loop(5), ind_loop(6)
           do  j = ind_loop(3), ind_loop(4)
            do  i = ind_loop(1), ind_loop(2)    
                              
              l = inddm(i,j,k)                            
        

               lstk  =  (i+1 - ind_loop(1))
     &                +(j-ind_loop(3))*nistk
     &             +(k-ind_loop(5))*nistk*nistk2
     &         + nzone*nistk*nistk2*nistk3


             drodm(l,ne) =  stock(lstk,ne)

                                
           enddo

         enddo
        enddo
       enddo

     
      end if

    


          
       end
