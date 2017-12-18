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

C==============================================================================
C Extrapole (recopie) les valeur de l'interieur du domaine
C==============================================================================
      SUBROUTINE extrapallbvaluesf(size, nfld, x, ni, nj, nk, 
     &     border)

      IMPLICIT NONE
#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E  size
      INTEGER_E  nfld
      INTEGER_E  ni, nj, nk
      INTEGER_E  border

C_OUT
      REAL_E     x(0:size-1, 1:nfld)

C_LOCAL
      INTEGER_E i, j, k, n, ind, ind2
      INTEGER_E ninj, iref, jref, kref
C==============================================================================

      ninj = ni*nj

      DO n = 1, nfld

C     I interfaces
         DO k = 1, nk
            iref = border+1
            DO j = 1, nj
               DO i = 1, border
                  ind = i-1 + (j-1)*ni + (k-1)*ninj
                  ind2 = iref-1 + (j-1)*ni + (k-1)*ninj
                  x(ind, n) = x(ind2, n)
               ENDDO
            ENDDO
            iref = ni-border
            DO j = 1, nj
               DO i = ni-border+1, ni
                  ind = i-1 + (j-1)*ni + (k-1)*ninj
                  ind2 = iref-1 + (j-1)*ni + (k-1)*ninj
                  x(ind, n) = x(ind2, n)
               ENDDO
            ENDDO
         ENDDO

C     J interfaces
         DO k = 1, nk
            jref = border+1
            DO j = 1, border
               DO i = 1, ni
                  ind = i-1 + (j-1)*ni + (k-1)*ninj
                  ind2 = i-1 + (jref-1)*ni + (k-1)*ninj
                  x(ind, n) = x(ind2, n)
               ENDDO
            ENDDO
            jref = nj-border
            DO j = nj-border+1, nj
               DO i = 1, ni
                  ind = i-1 + (j-1)*ni + (k-1)*ninj
                  ind2 = i-1 + (jref-1)*ni + (k-1)*ninj
                  x(ind, n) = x(ind2, n)
               ENDDO
            ENDDO
         ENDDO
         
C     K interfaces
         kref = border+1
         DO k = 1, border
            DO j = 1, nj
               DO i = 1, ni
                  ind = i-1 + (j-1)*ni + (k-1)*ninj
                  ind2 = i-1 + (j-1)*ni + (kref-1)*ninj
                  x(ind, n) = x(ind2, n)
               ENDDO
            ENDDO
         ENDDO
         kref = nk-border
         DO k = nk-border+1, nk
            DO j = 1, nj
               DO i = 1, ni
                  ind = i-1 + (j-1)*ni + (k-1)*ninj
                  ind2 = i-1 + (j-1)*ni + (kref-1)*ninj
                  x(ind, n) = x(ind2, n)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      END
