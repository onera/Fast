        E_Int* ind_fen = param_int + pt_bc + BC_FEN;
        E_Int  inc_bc[3];
        if ( idir <= 2 ) 
           {
            inc_bc[0] = ind_fen[3] - ind_fen[2] + 1; // nombre element de la fenetre dans la direction J
            inc_bc[1] = ind_fen[2]; // debut indice j
            inc_bc[2] = ind_fen[4]; // debut indice k
           }  
        else if ( idir <= 4 )
           { 
            inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
            inc_bc[1] = ind_fen[0]; // debut indice i
            inc_bc[2] = ind_fen[4]; // debut indice k
           }  
        else{           
            inc_bc[0] = ind_fen[1] - ind_fen[0] + 1; // nombre element de la fenetre dans la direction I
            inc_bc[1] = ind_fen[0]; // debut indice i
            inc_bc[2] = ind_fen[2]; // debut indice j
           }
