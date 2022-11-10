c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 39 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine extract_res(ndom, param_int, param_real, 
     &                       ind_loop,
     &                       drodm, vol, ro_res)
c***********************************************************************
c_U   USER : DANDOIS
c
c     ACT
c_A    Extraction du residu explicite
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, ind_loop(6), param_int(0:*)
 
      REAL_E   drodm( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E     vol( param_int(NDIMDX_MTR) )
      REAL_E  ro_res( param_int(NDIMDX) , param_int(NEQ) )

      REAL_E param_real(0:*)
 
C Var loc
      INTEGER_E  l,i,j,k,lij,lt, ltij,lvo
      REAL_E volinv

#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

      IF(param_int(EXTRACT_RES).eq.1) THEN

        IF(param_int(NEQ).eq.6) THEN

         if(param_int(ITYPZONE).ne.2) then !domaine 3d general, domaine 3d k homogene, 2d
#include  "FastS/Compute/loop_begin.for"
            volinv = 1./vol(lvo)

            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
            ro_res(l,6)= volinv*drodm(l,6) 
#include  "FastS/Compute/loop_end.for"

         else
            volinv = 1./vol(1)
#include  "FastS/Compute/loop_begin.for"
            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
            ro_res(l,6)= volinv*drodm(l,6) 
#include  "FastS/Compute/loop_end.for"

         endif
        ELSE !neq=5

         if (param_int(ITYPZONE).ne.2) then
#include  "FastS/Compute/loop_begin.for"
            volinv = 1./vol(lvo)

            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
#include  "FastS/Compute/loop_end.for"
         else
            volinv = 1./vol(1)
#include  "FastS/Compute/loop_begin.for"
            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
#include  "FastS/Compute/loop_end.for"
         endif
        ENDIF!neq

      ELSE  !itype=2)

        volinv = 1./param_real(DTC)

        IF(param_int(NEQ).eq.6) THEN

         if(param_int(ITYPZONE).ne.2) then !domaine 3d general, domaine 3d k homogene, 2d
#include  "FastS/Compute/loop_begin.for"
            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
            ro_res(l,6)= volinv*drodm(l,6) 
#include  "FastS/Compute/loop_end.for"

         else
            volinv = 1./vol(1)
#include  "FastS/Compute/loop_begin.for"
            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
            ro_res(l,6)= volinv*drodm(l,6) 
#include  "FastS/Compute/loop_end.for"

         endif
        ELSE !neq=5

         if (param_int(ITYPZONE).ne.2) then
#include  "FastS/Compute/loop_begin.for"
            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
#include  "FastS/Compute/loop_end.for"
         else
#include  "FastS/Compute/loop_begin.for"
            ro_res(l,1)= volinv*drodm(l,1) 
            ro_res(l,2)= volinv*drodm(l,2) 
            ro_res(l,3)= volinv*drodm(l,3) 
            ro_res(l,4)= volinv*drodm(l,4) 
            ro_res(l,5)= volinv*drodm(l,5)
#include  "FastS/Compute/loop_end.for"
         endif
        ENDIF!neq
      
      ENDIF
      end
