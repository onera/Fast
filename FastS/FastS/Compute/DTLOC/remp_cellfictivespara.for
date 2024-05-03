c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine remp_cellfictivespara(param_intDnr,param_intRcv,ropDnr, 
     &ropRcv,listDnr,listRcv,nbpts,pt_deb,pt_fin)
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

      INTEGER_E param_intDnr(0:*)
      INTEGER_E param_intRcv(0:*)

      INTEGER_E nbpts,pt_deb,pt_fin

      INTEGER_E listDnr(0:nbpts-1)
      INTEGER_E listRcv(0:nbpts-1)

      REAL_E ropDnr( param_intDnr(NDIMDX),param_intDnr(NEQ))
      REAL_E ropRcv( param_intRcv(NDIMDX),param_intRcv(NEQ))
  
 
C Var local
      INTEGER_E ne,neq,i,indRcv,indDnr


      neq=param_intDnr(NEQ)

      do ne=1,neq
         do i=pt_deb,pt_fin
         
            indRcv = listRcv(i)
            indDnr = listDnr(i)

            ropRcv(indRcv+1,ne)=ropDnr(indDnr+1,ne)

            !print*, indRcv+1 
           

            end do
          end do



       end
