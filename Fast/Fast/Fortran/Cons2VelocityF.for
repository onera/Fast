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
C Calcule la vitesse a partir des variables conservatives
C =============================================================================
      SUBROUTINE cons2velocity(cellt, ro, rou, rov, row, 
     &                         u, v, w)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN    
      INTEGER_E      cellt ! Total Cell Number
      REAL_E         ro(0:cellt-1) ! Density
      REAL_E         rou(0:cellt-1) ! MomentumX
      REAL_E         rov(0:cellt-1) ! MomentumY
      REAL_E         row(0:cellt-1) ! MomentumZ
C_OUT
      REAL_E         u(0:cellt-1) ! U
      REAL_E         v(0:cellt-1) ! V
      REAL_E         w(0:cellt-1) ! W

C_LOCAL
      INTEGER_E      i
      REAL_E         roi
C==============================================================================
!$OMP PARALLEL PRIVATE(i,roi)
!$OMP DO
      DO i = 0, cellt-1
         roi = ONE / ro(i)
         u(i) = rou(i) * roi
         v(i) = rov(i) * roi
         w(i) = row(i) * roi
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
