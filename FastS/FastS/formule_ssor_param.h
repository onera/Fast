C    adresse point courant pour tableau de la taille d'un domaine 
      integer*4 indssor, i_5,j_5,k_5,i_size,j_size

      indssor(i_5,j_5,k_5,i_size,j_size) = 1
     &      + (i_5-ind_loop_sdm(1)+param_int(NIJK+3))
     &      + (j_5-ind_loop_sdm(3)+param_int(NIJK+3))*i_size
     &      + (k_5-ind_loop_sdm(5)+param_int(NIJK+4))*i_size*j_size
