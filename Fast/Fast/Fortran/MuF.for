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
C Calcule la viscosite selon la loi de Sutherland
C =============================================================================
      SUBROUTINE viscosity(cellt, betas, Cs, temp, mu)

      IMPLICIT NONE

C==============================================================================
C_IN    
      INTEGER_E      cellt ! Total Cell Number
      REAL_E         betas      ! 1st Sutherland's Law constant 
      REAL_E         Cs         ! 2nd Sutherland's Law constant 
      REAL_E         temp(0:cellt-1) ! Temperature Field
C_OUT
      REAL_E         mu(0:cellt-1) ! Molecular Viscosity Coefficient

C_LOCAL
      INTEGER_E      i
      REAL_E         t     ! Local Temperature Value
C==============================================================================
!$OMP PARALLEL PRIVATE(i,t)
!$OMP DO
      DO i = 0, cellt-1
         t = temp(i)
         mu(i) = (betas * SQRT(t)*t) / (t + Cs)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
