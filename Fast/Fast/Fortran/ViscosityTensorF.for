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
      SUBROUTINE visctensor(
     &     cellt, mu,
     &     gradxu, gradyu, gradzu,
     &     gradxv, gradyv, gradzv,
     &     gradxw, gradyw, gradzw,
     &     tensxx, tensxy, tensxz,
     &     tensyy, tensyz, tenszz)
      IMPLICIT NONE
 
#include "Def/DefFortranConst.h"
 
C==============================================================================
C_IN    
      INTEGER_E  cellt ! Total Cell Number
      REAL_E     mu(0:cellt-1) ! "Generic" Viscosity Coefficient
      REAL_E     gradxu(0:cellt-1) ! x-component of the gradient of ux
      REAL_E     gradxv(0:cellt-1) ! x-component of the gradient of uy
      REAL_E     gradxw(0:cellt-1) ! x-component of the gradient of uz
      REAL_E     gradyu(0:cellt-1) ! y-component of the gradient of ux
      REAL_E     gradyv(0:cellt-1) ! y-component of the gradient of uy
      REAL_E     gradyw(0:cellt-1) ! y-component of the gradient of uz
      REAL_E     gradzu(0:cellt-1) ! z-component of the gradient of ux
      REAL_E     gradzv(0:cellt-1) ! z-component of the gradient of uy
      REAL_E     gradzw(0:cellt-1) ! z-component of the gradient of uz
C_OUT
      REAL_E     tensxx(0:cellt-1) ! Tensor - xx-Component
      REAL_E     tensxy(0:cellt-1) ! Tensor - xy-Component
      REAL_E     tensxz(0:cellt-1) ! Tensor - xz-Component
      REAL_E     tensyy(0:cellt-1) ! Tensor - yy-Component
      REAL_E     tensyz(0:cellt-1) ! Tensor - yz-Component
      REAL_E     tenszz(0:cellt-1) ! Tensor - zz-Component
C
C_LOCAL
      INTEGER_E  l              ! Loop cell index
      REAL_E     ds3            ! Constant (2./3.)
      REAL_E     mul, muds3
C------------------------------------------------------------------------------
      ds3 = TWO / THREE
!$OMP PARALLEL PRIVATE(l,mul,muds3)
!$OMP DO
      DO l = 0, cellt-1
         mul = mu(l)
         muds3 = mul*ds3
         tensxx(l) = muds3*(TWO*gradxu(l) - gradyv(l) - gradzw(l))
         tensyy(l) = muds3*(TWO*gradyv(l) - gradxu(l) - gradzw(l))
         tenszz(l) = muds3*(TWO*gradzw(l) - gradxu(l) - gradyv(l)) 
       
         tensxy(l) = mul*(gradyu(l) + gradxv(l))
         tensxz(l) = mul*(gradzu(l) + gradxw(l))
         tensyz(l) = mul*(gradzv(l) + gradyw(l)) 
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
