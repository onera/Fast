C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

C =============================================================================
C Calcul la pression aux centres en fonction des grandeurs conservatives
C =============================================================================
       SUBROUTINE pressure(cellt,
     &                     ro, rou, rov, row, roE,
     &                     gam,   
     &                     p)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E           cellt ! Total Cell Number 
      REAL_E              ro(0:cellt-1)
      REAL_E              rou(0:cellt-1)
      REAL_E              rov(0:cellt-1)
      REAL_E              row(0:cellt-1)
      REAL_E              roE(0:cellt-1)
      REAL_E              gam

C_OUT
      REAL_E              p(0:cellt-1) ! Pressure Field 

C_LOCAL
      INTEGER_E           l
      REAL_E              rovax2, rovay2, rovaz2
      REAL_E              roi, eint, gam1
c------------------------------------------------------------------------------

      gam1 = gam-ONE

!$OMP PARALLEL PRIVATE(l,roi,rovax2,rovay2,rovaz2,eint)
!$OMP DO 
      DO l = 0, cellt-1 
         roi = ONE / ro(l)
         rovax2 =  rou(l) * rou(l)   
         rovay2 =  rov(l) * rov(l)  
         rovaz2 =  row(l) * row(l)
         eint  =  roE(l)

         p(l) = gam1 *
     &        (eint-ONE_HALF*roi*(rovax2+rovay2+rovaz2))
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
