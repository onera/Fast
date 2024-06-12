C    adresse interface pour tableau vitesse entrainement
      INTEGER_E indven, i_4,j_4,k_4

      indven(i_4,j_4,k_4) =  1
     &           + (i_4+param_int(NIJK_VENT+3)-1)*param_int(NIJK_VENT)
     &           + (j_4+param_int(NIJK_VENT+3)-1)*param_int(NIJK_VENT+1)
     &           + (k_4+param_int(NIJK_VENT+4)-1)*param_int(NIJK_VENT+2)
