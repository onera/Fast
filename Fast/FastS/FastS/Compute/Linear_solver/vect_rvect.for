c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine vect_rvect(param_int, ind_loop, 
     &                      vect1, vect2, value,  normL2)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    selection preconditionneur a droite pour Krylov
c
c     VAL
c
c     I/O
c_/    vect1,2  value :  IN
c_/    normL2         : OUT
c
c***********************************************************************²
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ind_loop(6), param_int(0:*)

      REAL_E vect1( param_int(NDIMDX), param_int(NEQ) )
      REAL_E vect2( param_int(NDIMDX), param_int(NEQ) )
      REAL_E value, normL2

!var loc
      INTEGER_E k, j, lij, l, ltij, lt, lvo
      REAL_E tmp1, tmp2, tmp3, tmp4,tmp5,tmp6

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"

            tmp1 = vect1(l,1) - value*vect2(l,1)
            tmp2 = vect1(l,2) - value*vect2(l,2)
            tmp3 = vect1(l,3) - value*vect2(l,3)
            tmp5 = vect1(l,5) - value*vect2(l,5)
               
            vect1(l, 1) = tmp1
            vect1(l, 2) = tmp2
            vect1(l, 3) = tmp3
            vect1(l, 4) = 0.
            vect1(l, 5) = tmp5
               
            normL2 = normL2 + tmp1*tmp1 
     &                      + tmp2*tmp2 
     &                      + tmp3*tmp3 
     &                      + tmp5*tmp5
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
            tmp1 = vect1(l,1) - value*vect2(l,1)
            tmp2 = vect1(l,2) - value*vect2(l,2)
            tmp3 = vect1(l,3) - value*vect2(l,3)
            tmp5 = vect1(l,5) - value*vect2(l,5)
            tmp6 = vect1(l,6) - value*vect2(l,6)
               
            vect1(l, 1) = tmp1
            vect1(l, 2) = tmp2
            vect1(l, 3) = tmp3
            vect1(l, 4) = 0.
            vect1(l, 5) = tmp5
            vect1(l, 6) = tmp6
               
            normL2 = normL2 + tmp1*tmp1 
     &                      + tmp2*tmp2 
     &                      + tmp3*tmp3 
     &                      + tmp5*tmp5
     &                      + tmp6*tmp6
#include   "FastS/Compute/loop_end.for"
        endif

      !!
      !!zone 3d
      !!
      ELSE

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
            tmp1 = vect1(l,1) - value*vect2(l,1)
            tmp2 = vect1(l,2) - value*vect2(l,2)
            tmp3 = vect1(l,3) - value*vect2(l,3)
            tmp4 = vect1(l,4) - value*vect2(l,4)
            tmp5 = vect1(l,5) - value*vect2(l,5)
               
            vect1(l, 1) = tmp1
            vect1(l, 2) = tmp2
            vect1(l, 3) = tmp3
            vect1(l, 4) = tmp4
            vect1(l, 5) = tmp5
               
            normL2 = normL2 + tmp1*tmp1 
     &                      + tmp2*tmp2 
     &                      + tmp3*tmp3 
     &                      + tmp4*tmp4 
     &                      + tmp5*tmp5
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
            tmp1 = vect1(l,1) - value*vect2(l,1)
            tmp2 = vect1(l,2) - value*vect2(l,2)
            tmp3 = vect1(l,3) - value*vect2(l,3)
            tmp4 = vect1(l,4) - value*vect2(l,4)
            tmp5 = vect1(l,5) - value*vect2(l,5)
            tmp6 = vect1(l,6) - value*vect2(l,6)
               
            vect1(l, 1) = tmp1
            vect1(l, 2) = tmp2
            vect1(l, 3) = tmp3
            vect1(l, 4) = tmp4
            vect1(l, 5) = tmp5
            vect1(l, 6) = tmp6
               
            normL2 = normL2 + tmp1*tmp1 
     &                      + tmp2*tmp2 
     &                      + tmp3*tmp3 
     &                      + tmp4*tmp4 
     &                      + tmp5*tmp5
     &                      + tmp6*tmp6
#include   "FastS/Compute/loop_end.for"
        endif
      ENDIF

      end

