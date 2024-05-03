c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine ibcsource(ndom, param_int, ind_rhs,
     &                     rop, cellN, coe, drodm)
c***********************************************************************
c_P                          O N E R A
c
c
c     ACT
c_A    Mise a jour de drodm prise en compte IBC ordre zero a la funk
c
c     INP
c_I    ndom   : drodm   (rhs)
c_I    ipara  : cellN_IBC  (mask pour ibc)
c_I    rop    :  le champ ro, u,T
c     OUT
c
c     I/O
c_/    drodm  : increment des variables conservatives
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, param_int(0:*), ind_rhs(6)

      REAL_E drodm( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   coe( param_int(NDIMDX) * param_int(NEQ_COE) )
      REAL_E cellN( param_int(NDIMDX) )
      
 
C Var loc
      INTEGER_E incmax,l, i,j,k,ne,lij,ltij,lt,b,ind,vg, lvo,
     & ind_loop(6), v1,v2,v3,v4,v5,v6
      REAL_E ratio,coefH,xmut(1), rho,vu,vv,vw,q,ampli,dt !!ajout pour feinter option de vecto 
      
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
      
      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v6 = 5*param_int(NDIMDX)

      !prise en compte du masquage sous-domaine eventuel
       ind_loop(1) =  max( ind_rhs(1), param_int(IBC +1) )
       ind_loop(2) =  min( ind_rhs(2), param_int(IBC +2) )
       ind_loop(3) =  max( ind_rhs(3), param_int(IBC +3) )
       ind_loop(4) =  min( ind_rhs(4), param_int(IBC +4) )
       ind_loop(5) =  max( ind_rhs(5), param_int(IBC +5) )
       ind_loop(6) =  min( ind_rhs(6), param_int(IBC +6) )

      IF(param_int(NEQ).eq.5) THEN

        do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
#include    "FastS/Compute/loopI_begin.for"

             rho  = rop(l + v1)
             vu   = rop(l + v2)
             vv   = rop(l + v3)
             vw   = rop(l + v4)
             q    = 0.5*(vu*vu+vv*vv+vw*vw)            

! protection pour interaction CellN et CellNIBC
             dt = max(1e-30,coe(l +V1))
             
!ampli= cellN(l)/coe(l +V1)*rho
             ampli= cellN(l)/dt*rho
             
!     drodm(l +v1) = drodm(l +v1)*xcomp
             drodm(l +v2) = drodm(l +v2)+ vu*ampli
             drodm(l +v3) = drodm(l +v3)+ vv*ampli
             drodm(l +v4) = drodm(l +v4)+ vw*ampli
             drodm(l +v5) = drodm(l +v5)+  q*ampli
             
            enddo
          enddo
        enddo

      ELSE

        do  k = ind_loop(5), ind_loop(6)
          do  j = ind_loop(3), ind_loop(4)
#include    "FastS/Compute/loopI_begin.for"

             rho  = rop(l + v1)
             vu   = rop(l + v2)
             vv   = rop(l + v3)
             vw   = rop(l + v4)
             q    = 0.5*(vu*vu+vv*vv+vw*vw)

!     protection pour interaction CellN et CellNIBC
             dt = max(1e-30,coe(l +V1))

!ampli= cellN(l)/coe(l +V1)*rho
             ampli= cellN(l)/dt*rho
             
!drodm(l +v1) = drodm(l +v1)*xcomp
             drodm(l +v2) = drodm(l +v2)+ vu*ampli
             drodm(l +v3) = drodm(l +v3)+ vv*ampli
             drodm(l +v4) = drodm(l +v4)+ vw*ampli
             drodm(l +v5) = drodm(l +v5)+  q*ampli
             drodm(l +v6) = drodm(l +v6)+ rop(l + v6)*ampli

            enddo
          enddo
        enddo

      ENDIF

      end
