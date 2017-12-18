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
      SUBROUTINE densfluxdiff(
     &     cellt,
     &     tgenxx, tgenxy, tgenxz,
     &     tgenyy, tgenyz, tgenzz,
     &     vax, vay, vaz, 
     &     qgenx,  qgeny,  qgenz,
     &     fmomxdx, fmomxdy, fmomxdz,
     &     fmomydx, fmomydy, fmomydz,
     &     fmomzdx, fmomzdy, fmomzdz,
     &     froedx,  froedy,  froedz)

      IMPLICIT NONE 
 
C==============================================================================
C_IN    
      INTEGER_E  cellt ! Total Cell Number
      REAL_E     tgenxx(0:cellt-1) ! Tensor - xx-Component
      REAL_E     tgenxy(0:cellt-1) ! Tensor - xy-Component
      REAL_E     tgenxz(0:cellt-1) ! Tensor - xz-Component
      REAL_E     tgenyy(0:cellt-1) ! Tensor - yy-Component
      REAL_E     tgenyz(0:cellt-1) ! Tensor - yz-Component
      REAL_E     tgenzz(0:cellt-1) ! Tensor - zz-Component
      REAL_E     vax(0:cellt-1) ! Velocity - x-Component
      REAL_E     vay(0:cellt-1) ! Velocity - y-Component
      REAL_E     vaz(0:cellt-1) ! Velocity - z-Component
      REAL_E     qgenx(0:cellt-1) ! Heat Flux Vector - x-Component
      REAL_E     qgeny(0:cellt-1) ! Heat Flux Vector - y-Component
      REAL_E     qgenz(0:cellt-1) ! Heat Flux Vector - z-Component
C_OUT
      REAL_E     fmomxdx(0:cellt-1) ! Diffusive Flux Density for MomentumX - x-Component
      REAL_E     fmomxdy(0:cellt-1) ! Diffusive Flux Density for MomentumX - y-Component
      REAL_E     fmomxdz(0:cellt-1) ! Diffusive Flux Density for MomentumX - z-Component
      REAL_E     fmomydx(0:cellt-1) ! Diffusive Flux Density for MomentumY - x-Component
      REAL_E     fmomydy(0:cellt-1) ! Diffusive Flux Density for MomentumY - y-Component
      REAL_E     fmomydz(0:cellt-1) ! Diffusive Flux Density for MomentumY - z-Component
      REAL_E     fmomzdx(0:cellt-1) ! Diffusive Flux Density for MomentumZ - x-Component
      REAL_E     fmomzdy(0:cellt-1) ! Diffusive Flux Density for MomentumZ - y-Component
      REAL_E     fmomzdz(0:cellt-1) ! Diffusive Flux Density for MomentumZ - z-Component
      REAL_E     froedx(0:cellt-1) ! Diffusive Flux Density for Energy - x-Component
      REAL_E     froedy(0:cellt-1) ! Diffusive Flux Density for Energy - y-Component
      REAL_E     froedz(0:cellt-1) ! Diffusive Flux Density for Energy - z-Component
C_LOCAL
      INTEGER_E  i
      REAL_E txx, txy, txz, tyy, tyz, tzz, vx, vy, vz
C -----------------------------------------------------------------------------
!$OMP PARALLEL PRIVATE(i,txx,tyy,tzz,txy,txz,tyz,vx,vy,vz)
!$OMP DO
      DO i = 0, cellt-1
         txx = tgenxx(i)
         tyy = tgenyy(i)
         tzz = tgenzz(i)
         txy = tgenxy(i)
         txz = tgenxz(i)
         tyz = tgenyz(i)
         vx = vax(i)
         vy = vay(i)
         vz = vaz(i)

         fmomxdx(i) = -txx
         fmomxdy(i) = -txy
         fmomxdz(i) = -txz

         fmomydx(i) = -txy
         fmomydy(i) = -tyy
         fmomydz(i) = -tyz 
         
         fmomzdx(i) = -txz
         fmomzdy(i) = -tyz
         fmomzdz(i) = -tzz
         
         froedx(i) = -txx*vx-txy*vy-txz*vz+qgenx(i) 
         froedy(i) = -txy*vx-tyy*vy-tyz*vz+qgeny(i) 
         froedz(i) = -txz*vx-tyz*vy-tzz*vz+qgenz(i)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
