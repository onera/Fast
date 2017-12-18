C    adresse interface pour tableau metric
      integer*4 indmtr, i_3,j_3,k_3

      indmtr(i_3,j_3,k_3) =  1
     &                   + (i_3+param_int(NIJK_MTR+3)-1)*param_int(NIJK_MTR)
     &                   + (j_3+param_int(NIJK_MTR+3)-1)*param_int(NIJK_MTR+1)
     &                   + (k_3+param_int(NIJK_MTR+4)-1)*param_int(NIJK_MTR+2)
