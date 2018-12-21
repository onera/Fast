C  
C    Copyright 2013-2019 Onera.
C
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
C Local time step computation, convective part
C =============================================================================
      SUBROUTINE cptimestepconv(cellNbTot, cellDim,
     &     ro, rou, rov, row, roE,
     &     dtSteady,
     &     gamma, cfl)
      IMPLICIT NONE
# include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E           cellNbTot
      REAL_E              ro(0:cellNbTot-1)
      REAL_E              rou(0:cellNbTot-1)
      REAL_E              rov(0:cellNbTot-1)
      REAL_E              row(0:cellNbTot-1)
      REAL_E              roE(0:cellNbTot-1)
      REAL_E              gamma, cfl, cellDim

C_OUT
      REAL_E              dtSteady(0:cellNbTot-1)
C_LOCAL
      INTEGER_E           i
      REAL_E              rhoi, velo2, velo, gamma1
      REAL_E              ener, sound, ccc
      REAL_E              g2
c------------------------------------------------------------------------------

      ccc = cfl * cellDim
      gamma1 = gamma-ONE
      g2  = gamma * gamma1

!$OMP PARALLEL PRIVATE(i,rhoi,velo2,velo,ener,sound)
!$OMP DO 
      DO i = 0, cellNbTot-1
        rhoi   = ONE / ro(i)
        velo2 = rou(i) * rou(i) +
     &       rov(i) * rov(i) +
     &       row(i) * row(i)
        velo  = SQRT(velo2) * rhoi
        ener  = roE(i) - ONE_HALF * velo2 * rhoi
        sound = SQRT(g2 * ener * rhoi)

        dtSteady(i) = ccc / (velo + sound)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
