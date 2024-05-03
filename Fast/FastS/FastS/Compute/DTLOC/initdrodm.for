c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine initdrodm(param_int,ind_loop,drodm,drodm2)
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

      INTEGER_E ind_loop(6), param_int(0:*)

      !REAL_E drodm(param_int(NDIMDX) ,param_int(NEQ))
      !REAL_E drodm2(param_int(NDIMDX) ,param_int(NEQ))

      REAL_E drodm(param_int(NDIMDX) * param_int(NEQ) )
      REAL_E drodm2(param_int(NDIMDX) * param_int(NEQ))
    
C Var local

      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,neq,n,lstk,vg,lvo
      REAL_E ratio,coefH,xmut(1),rop(1)
      REAL_E c1,a

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
      
      neq=param_int(NEQ)

 
      !print*, ind_loop(1), ind_loop(2), ind_loop(3), ind_loop(4)

c$$$      do  ne=1,neq
c$$$       do  k = ind_loop(5), ind_loop(6)
c$$$        do  j = ind_loop(3), ind_loop(4)
c$$$         do  i = ind_loop(1), ind_loop(2)                              
c$$$               
c$$$              l  = inddm(i,j,k)
c$$$
c$$$              drodm2(l,ne) = drodm(l,ne)
c$$$              
c$$$              !if (i.le.ind_loop(1)+6 .and. j==50 .and. ne==1) then
c$$$              !   print*,"drodm2= ",drodm2(l,ne),k 
c$$$              !end if
c$$$              
c$$$
c$$$           enddo
c$$$                          
c$$$         enddo
c$$$         enddo
c$$$         enddo


         do  ne= 1, neq
          vg = param_int(NDIMDX)*(ne-1)
         do  k = ind_loop(5), ind_loop(6)
         do  j = ind_loop(3), ind_loop(4)
#include    "FastS/Compute/loopI_begin.for"

             
             drodm2(l+vg) = drodm(l+vg)
              
              !if (i.le.ind_loop(1)+6 .and. j==50 .and. ne==1) then
              !   print*,"drodm2= ",drodm2(l,ne),k 
              !end if
              

           enddo
                          
         enddo
         enddo
         enddo




          
         
       end
