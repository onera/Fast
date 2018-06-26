***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine dp_dw_vect(param_int, param_real, ind_loop, rop,
     &                      vectin, vectout)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    Jacobienne dp_dw: 
c      p               : (density,Velocity,Temperature, nutilde)
c      w               : (density,Momentum,TotalEnergy, Densitynutilde)
c
c     VAL
c
c     I/O
c_/    vectin    : vecteur kryloc IN
c_/    vectout   : vecteur kryloc OUT
c
c***********************************************************************      
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E param_int(0:*), ind_loop(6)
      REAL_E param_real(0:*)
      REAL_E     rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E  vectin( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E vectout( param_int(NDIMDX) , param_int(NEQ) )

      INTEGER_E k, j, lij, l, lt,ltij,lmtr,lvo
      REAL_E cvinv, roinv, ke, 
     &     dP2_dW1, dP2_dW2, dP3_dW1, dP3_dW3, dP4_dW1, dP4_dW4, 
     &     dP5_dW1, dP5_dW2, dP5_dW3, dP5_dW4, dP5_dW5,
     &     dP6_dW1, dP6_dW6

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      cvinv = 1./param_real(CVINF)


      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"

               ke   = 0.5 * ( rop(l,2)*rop(l,2) +rop(l,3)*rop(l,3))

               roinv   = 1./rop(l,1)

               dP2_dW1 = - rop(l,2)*roinv
               dP2_dW2 = roinv

               dP3_dW1 = - rop(l,3)*roinv
               dP3_dW3 = roinv

               dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
               dP5_dW2 = dP2_dW1*cvinv
               dP5_dW3 = dP3_dW1*cvinv
               dP5_dW5 = roinv*cvinv 

               vectout(l,1) = vectin(l,1)
               vectout(l,2) = dP2_dW1 *vectin(l,1) + dP2_dW2*vectin(l,2)
               vectout(l,3) = dP3_dW1 *vectin(l,1) + dP3_dW3*vectin(l,3)
               vectout(l,4) = 0.
               vectout(l,5) = dP5_dW1 *vectin(l,1) + dP5_dW2*vectin(l,2)
     &                       +dP5_dW3 *vectin(l,3) + dP5_dW5*vectin(l,5)
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"

               ke   = 0.5 * ( rop(l,2)*rop(l,2) +rop(l,3)*rop(l,3))

               roinv   = 1./rop(l,1)

               dP2_dW2 = roinv
               dP2_dW1 = - rop(l,2)*roinv

               dP3_dW1 = - rop(l,3)*roinv
               dP3_dW3 = roinv

               dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
               dP5_dW2 = dP2_dW1*cvinv
               dP5_dW3 = dP3_dW1*cvinv
               dP5_dW5 = roinv*cvinv 

               dP6_dW1 = - rop(l,6)*roinv
               dP6_dW6 = roinv

               vectout(l,1) = vectin(l,1)
               vectout(l,2) = dP2_dW1 *vectin(l,1) + dP2_dW2*vectin(l,2)
               vectout(l,3) = dP3_dW1 *vectin(l,1) + dP3_dW3*vectin(l,3)
               vectout(l,4) = 0.
               vectout(l,5) = dP5_dW1 *vectin(l,1) + dP5_dW2*vectin(l,2)
     &                       +dP5_dW3 *vectin(l,3) + dP5_dW5*vectin(l,5)
               vectout(l,6) = dP6_dW1 *vectin(l,1) + dP6_dW6*vectin(l,6)
#include   "FastS/Compute/loop_end.for"
        endif

      !!
      !!zone 3d
      !!
      ELSE


        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"

               ke   = 0.5*(  rop(l,2)*rop(l,2)
     &                     + rop(l,3)*rop(l,3)
     &                     + rop(l,4)*rop(l,4))

               roinv   = 1./rop(l,1)

               dP2_dW1 = - rop(l,2)*roinv
               dP2_dW2 = roinv

               dP3_dW1 = - rop(l,3)*roinv
               dP3_dW3 = roinv

               dP4_dW1 = - rop(l,4)*roinv
               dP4_dW4 = roinv

               dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
               dP5_dW2 = dP2_dW1*cvinv
               dP5_dW3 = dP3_dW1*cvinv
               dP5_dW4 = dP4_dW1*cvinv
               dP5_dW5 = roinv*cvinv 

               vectout(l,1) = vectin(l,1)
               vectout(l,2) = dP2_dW1*vectin(l,1) + dP2_dW2*vectin(l,2)
               vectout(l,3) = dP3_dW1*vectin(l,1) + dP3_dW3*vectin(l,3)
               vectout(l,4) = dP3_dW1*vectin(l,1) + dP4_dW4*vectin(l,4)
               vectout(l,5) = dP5_dW1*vectin(l,1) + dP5_dW2*vectin(l,2)
     &                       +dP5_dW3*vectin(l,3) + dP5_dW4*vectin(l,4)
     &                       +dP5_dW5*vectin(l,5)
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"

               ke   = 0.5*(  rop(l,2)*rop(l,2)
     &                     + rop(l,3)*rop(l,3)
     &                     + rop(l,4)*rop(l,4))

               roinv   = 1./rop(l,1)

               dP2_dW1 = - rop(l,2)*roinv
               dP2_dW2 = roinv

               dP3_dW1 = - rop(l,3)*roinv
               dP3_dW3 = roinv

               dP4_dW1 = - rop(l,4)*roinv
               dP4_dW4 = roinv

               dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
               dP5_dW2 = dP2_dW1*cvinv
               dP5_dW3 = dP3_dW1*cvinv
               dP5_dW4 = dP4_dW1*cvinv
               dP5_dW5 = roinv*cvinv 

               vectout(l,1) = vectin(l,1)
               vectout(l,2) = dP2_dW1*vectin(l,1) + dP2_dW2*vectin(l,2)
               vectout(l,3) = dP3_dW1*vectin(l,1) + dP3_dW3*vectin(l,3)
               vectout(l,4) = dP3_dW1*vectin(l,1) + dP4_dW4*vectin(l,4)
               vectout(l,5) = dP5_dW1*vectin(l,1) + dP5_dW2*vectin(l,2)
     &                       +dP5_dW3*vectin(l,3) + dP5_dW4*vectin(l,4)
     &                       +dP5_dW5*vectin(l,5)
               vectout(l,6) = dP6_dW1*vectin(l,1) + dP6_dW6*vectin(l,6)
#include   "FastS/Compute/loop_end.for"
        endif

      ENDIF

      end
