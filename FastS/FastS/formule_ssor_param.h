C    adresse point courant pour tableau de la taille d'un domaine 
      integer*4 indssor, i_5,j_5,k_5

      indssor(i_5,j_5,k_5) = 1  
     &     + (i_5+param_int(NIJK+3)-1)
     &     + (j_5+param_int(NIJK+3)-1)*param_int(NIJK)
     &     + (k_5+param_int(NIJK+4)-1)*param_int(NIJK)*param_int(NIJK+1)
