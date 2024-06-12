C    adresse point courant pour tableau de la taille d'un domaine 
      INTEGER_E indssor, i_5,j_5,k_5,ssor_i,ssor_j

      indssor(i_5,j_5,k_5) = 1
     &      + (i_5-ind_loop_sdm(1)+param_int(NIJK+3))
     &      + (j_5-ind_loop_sdm(3)+param_int(NIJK+3))*ssor_i
     &      + (k_5-ind_loop_sdm(5)+param_int(NIJK+4))*ssor_i*ssor_j
