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
C    
C =============================================================================
C Calcule la temperature a partir des variables conservatives
C =============================================================================
      SUBROUTINE temp(cellt, cv, ro, u, v, w, roE,
     &                t)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN    
      INTEGER_E      cellt ! Total Cell Number
      REAL_E         cv
      REAL_E         ro(0:cellt-1) ! Density
      REAL_E         u(0:cellt-1) ! MomentumX
      REAL_E         v(0:cellt-1) ! MomentumY
      REAL_E         w(0:cellt-1) ! MomentumZ
      REAL_E         roE(0:cellt-1) ! MomentumZ
C_OUT
      REAL_E         t(0:cellt-1) ! temperature

C_LOCAL
      INTEGER_E      i
      REAL_E         roi, cvi, u2, v2, w2
      
C==============================================================================
      cvi = ONE / cv
!$OMP PARALLEL PRIVATE(i,u2,v2,w2)
!$OMP DO
      DO i = 0, cellt-1
         u2 = u(i)
         u2 = u2*u2
         v2 = v(i)
         v2 = v2*v2
         w2 = w(i)
         w2 = w2*w2
         t(i) = cvi*(roE(i)/ro(i)-ONE_HALF*(u2+v2+w2))
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
