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
C Local time step computation, turbulent part
C =============================================================================
      SUBROUTINE cptimesteptur(cellNbTot, ro,
     &     cellDimension,
     &     dtSteady, timeStepConv,
     &     mu, muTurb,
     &     gamma, prandtl, pranTurb, cfl)
      IMPLICIT NONE
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E           cellNbTot
      REAL_E              gamma, prandtl, pranTurb, cfl
      REAL_E              ro(0:cellNbTot-1)
      REAL_E              cellDimension(1)
      REAL_E              timeStepConv(0:cellNbTot-1)
      REAL_E              mu(0:cellNbTot-1), muTurb(0:cellNbTot-1)

C_OUT
      REAL_E              dtSteady(0:cellNbTot-1)
C_LOCAL
      INTEGER_E           i
      REAL_E              rho, denom, dim, timeIncrV, ccc
      REAL_E              ptli, ptlti
c------------------------------------------------------------------------------
      dim = cellDimension(1)
      ccc = cfl * ONE_HALF * dim * dim

      ptli = gamma / prandtl
      ptlti = gamma / pranTurb

      DO i = 0, cellNbTot-1
        rho   = ro(i)
        denom =   mu(i)  * ptli
     &       + muTurb(i) * ptlti
       
        timeIncrV = rho / denom
        dtSteady(i) =  MIN(timeStepConv(i), timeIncrV)
      END DO

      END
