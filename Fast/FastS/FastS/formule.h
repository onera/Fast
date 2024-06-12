C    adresse point courant pour tableau de la taille d'un domaine 
      INTEGER_E inddm, i_1,j_1,k_1

      inddm(i_1,j_1,k_1) = 1  
     &                   + (i_1+nijk(4)-1)
     &                   + (j_1+nijk(4)-1)*nijk(1)
     &                   + (k_1+nijk(5)-1)*nijk(1)*nijk(2)
