C    adresse interface pour tableau metric
      INTEGER_E indmtr, i_3,j_3,k_3

      indmtr(i_3,j_3,k_3) =  1
     &                   + (i_3+nijk_mtr(4)-1)*nijk_mtr(1)
     &                   + (j_3+nijk_mtr(4)-1)*nijk_mtr(2)
     &                   + (k_3+nijk_mtr(5)-1)*nijk_mtr(3)
