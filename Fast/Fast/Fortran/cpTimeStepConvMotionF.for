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
C Local time step computation, convective part with motion
C ============================================================================
      SUBROUTINE cptimestepconvmotion(cellNbTot, 
     &     sx, sy, sz,
     &     cellDimension,
     &     ro, rou, rov, row, roE, 
     &     dtSteady,
     &     gamma, cfl)
      IMPLICIT NONE
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E           cellNbTot
      REAL_E              ro(0:cellNbTot-1)
      REAL_E              rou(0:cellNbTot-1)
      REAL_E              rov(0:cellNbTot-1)
      REAL_E              row(0:cellNbTot-1)
      REAL_E              roE(0:cellNbTot-1)
      REAL_E              gamma, cfl
      REAL_E              sx(0:cellNbTot-1)
      REAL_E              sy(0:cellNbTot-1)
      REAL_E              sz(0:cellNbTot-1)
      REAL_E              cellDimension(1)

C_OUT
      REAL_E              dtSteady(0:cellNbTot-1)
C_LOCAL
      INTEGER_E           i
      REAL_E              rhoi, velo2, veloX, veloY, veloZ
      REAL_E              ener, sound, ccc
      REAL_E              veloMotionX, veloMotionY, veloMotionZ
      REAL_E              veloMotion, g2, gamma1
c-----------------------------------------------------------------------------

      ccc = cfl * cellDimension(1)
      gamma1 = gamma-ONE
      g2  = gamma * gamma1

      DO i = 0, cellNbTot-1
        rhoi   = ONE / ro(i)
        veloX = rou(i) * rhoi
        veloY = rov(i) * rhoi
        veloZ = row(i) * rhoi

        veloMotionX = veloX - sx(i)
        veloMotionY = veloY - sy(i)
        veloMotionZ = veloZ - sz(i)
        
        velo2 = veloX*veloX + veloY*veloY + veloZ*veloZ

        veloMotion = SQRT(veloMotionX*veloMotionX +
     &                    veloMotionY*veloMotionY +
     &                    veloMotionZ*veloMotionZ)

        ener  = roE(i) * rhoi - ONE_HALF * velo2
        sound = SQRT(g2 * ener)

        dtSteady(i) = ccc / (veloMotion + sound)
      END DO

      END
