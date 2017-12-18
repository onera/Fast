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
C Calcule le coefficient de conductivite thermique
C =============================================================================
      SUBROUTINE heatcoef(cellt, cp, prandtl, mu,
     &                    kappa)

      IMPLICIT NONE

C==============================================================================
C_IN    
      INTEGER_E      cellt ! Total Cell Number
      REAL_E         mu(0:cellt-1) ! viscosity
      REAL_E         cp
      REAL_E         prandtl
C_OUT
      REAL_E         kappa(0:cellt-1) ! Heat conductivity coefficient

C_LOCAL
      REAL_E        CpSPr
      INTEGER_E     i
C==============================================================================
      CpSPr = cp/prandtl
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = 0, cellt-1
        kappa(i) = CpSPr * mu(i)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
