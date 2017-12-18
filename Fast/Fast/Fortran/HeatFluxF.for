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
C Calcule le flux de chaleur (loi de Fourier)
C =============================================================================
      SUBROUTINE heatflux(cellt, kappa, gradxtemp, gradytemp, gradztemp,
     &                    qx, qy, qz)

      IMPLICIT NONE

C==============================================================================
C_IN    
      INTEGER_E      cellt ! Total Cell Number
      REAL_E         gradxtemp(0:cellt-1) ! gradient de temperature (1e composante)
      REAL_E         gradytemp(0:cellt-1) ! gradient de temperature (2e composante)
      REAL_E         gradztemp(0:cellt-1) ! gradient de temperature (3e composante)

      REAL_E         kappa(0:cellt-1) ! Heat conductivity coefficient
C_OUT
      REAL_E         qx(0:cellt-1) ! 1st Heat flux component
      REAL_E         qy(0:cellt-1) ! 2nd Heat flux component
      REAL_E         qz(0:cellt-1) ! 3rd Heat flux component

C_LOCAL
      INTEGER_E     i
      REAL_E        c
C==============================================================================
!$OMP PARALLEL PRIVATE(i,c)
!$OMP DO
      DO i = 0, cellt-1
         c = -kappa(i)
         qx(i) = c * gradxtemp(i)
         qy(i) = c * gradytemp(i)
         qz(i) = c * gradztemp(i) 
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
