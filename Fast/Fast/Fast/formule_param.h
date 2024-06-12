C    adresse point courant pour tableau de la taille d'un domaine 
      INTEGER_E inddm, i_1,j_1,k_1

      inddm(i_1,j_1,k_1) = 1  
     &     + (i_1+param_int(NIJK+3)-1)
     &     + (j_1+param_int(NIJK+3)-1)*param_int(NIJK)
     &     + (k_1+param_int(NIJK+4)-1)*param_int(NIJK)*param_int(NIJK+1)
