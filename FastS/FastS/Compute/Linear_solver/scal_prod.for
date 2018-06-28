c***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine scal_prod(param_int, ind_loop, vect1, vect2, value)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    calcul produit scalaire vec1.vec2 stocké dans value
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

      INTEGER_E ind_loop(6), param_int(0:*)


      REAL_E vect1( param_int(NDIMDX), param_int(NEQ) )
      REAL_E vect2( param_int(NDIMDX), param_int(NEQ) )
      REAL_E value

      INTEGER_E k, j, lij, l, ltij, lt, lvo

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
               value = value + vect1(l,1) * vect2(l,1)
     &                       + vect1(l,2) * vect2(l,2)
     &                       + vect1(l,3) * vect2(l,3) 
     &                       + vect1(l,5) * vect2(l,5)
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
               value = value + vect1(l,1) * vect2(l,1)
     &                       + vect1(l,2) * vect2(l,2)
     &                       + vect1(l,3) * vect2(l,3) 
     &                       + vect1(l,5) * vect2(l,5)
     &                       + vect1(l,6) * vect2(l,6)
#include   "FastS/Compute/loop_end.for"
        endif

      !!
      !!zone 3d
      !!
      ELSE

        if(param_int(NEQ).eq.5) then
#include   "FastS/Compute/loop_begin.for"
               value = value + vect1(l,1) * vect2(l,1)
     &                       + vect1(l,2) * vect2(l,2)
     &                       + vect1(l,3) * vect2(l,3) 
     &                       + vect1(l,4) * vect2(l,4) 
     &                       + vect1(l,5) * vect2(l,5)
#include   "FastS/Compute/loop_end.for"
        else
#include   "FastS/Compute/loop_begin.for"
               value = value + vect1(l,1) * vect2(l,1)
     &                       + vect1(l,2) * vect2(l,2)
     &                       + vect1(l,3) * vect2(l,3) 
     &                       + vect1(l,4) * vect2(l,4) 
     &                       + vect1(l,5) * vect2(l,5)
     &                       + vect1(l,6) * vect2(l,6)
#include   "FastS/Compute/loop_end.for"
        endif
      ENDIF

      end
