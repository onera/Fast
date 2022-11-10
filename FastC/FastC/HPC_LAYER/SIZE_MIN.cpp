      //Calcul de la taille minimal 1D du bloc thread 
      if (param_int[nd][IFLOW]==4)      { lmin=1;} //lbm
      else if(param_int[nd][ITYPCP]==2) { lmin=4;} // ns explicit
      else {lmin =10;}
