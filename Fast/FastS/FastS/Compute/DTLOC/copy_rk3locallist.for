c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine copy_rk3locallist(param_int, rop, stock, listDnr,
     & nbpts)
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

      INTEGER_E param_int(0:*), nbpts
      INTEGER_E listDnr(nbpts)

      REAL_E rop( param_int(NDIMDX),param_int(NEQ))
      REAL_E stock(nbpts,param_int(NEQ))     
 
C Var local
      INTEGER_E i, indDnr, ne

        
         do  ne=1,param_int(NEQ)
             do  i = 1, nbpts
                              
                indDnr = listDnr(i)

                stock(i,ne) = rop(indDnr+1,ne)

             end do
          end do           

      
      end
