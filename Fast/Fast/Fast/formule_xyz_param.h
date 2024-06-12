C    adresse point courant pour tableau x ou de la taille d'un domaine 
      INTEGER_E indcg, i_2,j_2,k_2

      indcg(i_2,j_2,k_2) = 1  
     &     + (i_2+param_int(NIJK_XYZ+3)-1)
     &     + (j_2+param_int(NIJK_XYZ+3)-1)*param_int(NIJK_XYZ)
     &     + (k_2+param_int(NIJK_XYZ+4)-1)
     *       *param_int(NIJK_XYZ)*param_int(NIJK_XYZ+1)
