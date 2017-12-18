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
C Calcul la densite de flux convectifs aux centres
C =============================================================================
      SUBROUTINE denconvec(cellt, 
     &                     p,     
     &                     ro, rou, rov, row, roE,  
     &                     fcdx, fcdy, fcdz)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E  cellt
      REAL_E     p(0:cellt-1) ! Pressures
      REAL_E     ro(0:cellt-1)
      REAL_E     rou(0:cellt-1)
      REAL_E     rov(0:cellt-1)
      REAL_E     row(0:cellt-1)
      REAL_E     roE(0:cellt-1)
C_OUT
      REAL_E     fcdx(0:cellt-1,5) ! Convective Flux  X-Component 
      REAL_E     fcdy(0:cellt-1,5) ! Convective Flux  Y-Component
      REAL_E     fcdz(0:cellt-1,5) ! Convective Flux  Z-Component

C_LOCAL
      INTEGER_E  i     
      REAL_E     roi,enthal ! Inverse of Density 
      REAL_E     tu, tv, tw, pi
      REAL_E     trou, trov, trow
C==============================================================================

!$OMP PARALLEL PRIVATE(i,roi,pi,trou,trov,trow,enthal,tu,tv,tw)

!$OMP DO
      DO i = 0, cellt-1
         fcdx(i, 1)  = rou(i)
         fcdy(i, 1)  = rov(i)
         fcdz(i, 1)  = row(i)
         roi = ONE / ro(i)
         pi = p(i)
         trou = rou(i)
         trov = rov(i)
         trow = row(i)
         enthal =  roE(i) + pi
C     
         tu = trou * roi
         tv = trov * roi
         tw = trow * roi
C
         fcdx(i, 2) = trou * tu + pi
         fcdx(i, 3) = trov * tu 
         fcdx(i, 4) = trow * tu
         fcdy(i, 2) = fcdx(i, 3)
         fcdy(i, 3) = trov * tv + pi
         fcdy(i, 4) = trow * tv
         fcdz(i, 2) = fcdx(i, 4)
         fcdz(i, 3) = fcdy(i, 4)
         fcdz(i, 4) = trow * tw + pi 
C
         fcdx(i, 5) = enthal * tu
         fcdy(i, 5) = enthal * tv
         fcdz(i, 5) = enthal * tw
      END DO
!$OMP END DO
!$OMP END PARALLEL
      RETURN
      END
