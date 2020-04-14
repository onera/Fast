C    adresse point courant pour tableau x ou de la taille d'un domaine 
      integer*4 indcg, i_2,j_2,k_2

      indcg(i_2,j_2,k_2) =  1
     &                   + (i_2+nijk_xyz(4)-1)
     &                   + (j_2+nijk_xyz(4)-1)*nijk_xyz(1)
     &                   + (k_2+nijk_xyz(5)-1)*nijk_xyz(1)*nijk_xyz(2)
