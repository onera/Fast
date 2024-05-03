c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 59 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine tstb3_global(ndom,flagCellN, param_int, param_real,
     &                        ind_loop,
     &                        coe, vol, cellN)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Calcul du pas de temps par noeud en gaz parfait.
c
c     VAL
c_V       Valide pour calcul Euler/NS/NSturb 
c
c     OUT
c_O    coe(:,11)   : tableau des pas de temps = dtc /vol(ijk)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, flagCellN, param_int(0:*), ind_loop(6)

      REAL_E cellN( param_int(NDIMDX)     )
      REAL_E   vol( param_int(NDIMDX_MTR) ) 
      REAL_E   coe( param_int(NDIMDX) , param_int(NEQ_COE) )

      REAL_E param_real(0:*)
C Var loc
      INTEGER_E  l,i,j,k,lij,lt, ltij,lvo
      REAL_E volinv

#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

      If(flagCellN.eq.1) then  !zone chimere

         if(param_int(ITYPZONE).ne.2) then !domaine 3d general, domaine 3d k homogene, 2d

#include     "FastS/Compute/loop_begin.for" 
             coe(l,1)=param_real(DTC)/vol(lvo)*MIN(cellN(l),2.-cellN(l))
#include     "FastS/Compute/loop_end.for" 
           
         else  
           volinv = param_real(DTC)/vol(1)

#include     "FastS/Compute/loop_begin.for" 
              coe(l,1) = volinv*MIN(cellN(l),2.-cellN(l))
#include     "FastS/Compute/loop_end.for" 
         endif

      Else  ! non chimere

         if(param_int(ITYPZONE).ne.2) then !domaine 3d general, domaine 3d k homogene, 2d

#include     "FastS/Compute/loop_begin.for" 
               coe(l,1) = param_real(DTC)/vol(lvo)
#include     "FastS/Compute/loop_end.for" 
           
         else            
           volinv = param_real(DTC)/vol(1)

#include     "FastS/Compute/loop_begin.for" 
               coe(l,1) = volinv
#include     "FastS/Compute/loop_end.for" 
         endif

      Endif !chimere ou pas

      end
