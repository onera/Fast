                E_Int pos_tab = param_int_tc[debut_rac + 26];
		E_Int donorPts_[6]; 
		donorPts_[0] =  param_int_tc[debut_rac + 7];
		donorPts_[1] =  param_int_tc[debut_rac + 8] ;
		donorPts_[2] =  param_int_tc[debut_rac + 9];
		donorPts_[3] =  param_int_tc[debut_rac + 10];
		donorPts_[4] =  param_int_tc[debut_rac + 11];
		donorPts_[5] =  param_int_tc[debut_rac + 12];

                // pourquoi int et pas E_Int
		int dir = param_int_tc[debut_rac + 13];
		E_Int profondeur = param_int_tc[debut_rac + 20];

		donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] = donorPts_[2*abs(dir)-1-(dir+abs(dir))/(2*abs(dir))] -(dir/abs(dir))*(profondeur);

		E_Int taillefenetre;
		taillefenetre = (donorPts_[1] - donorPts_[0] + 1)*(donorPts_[3] - donorPts_[2] + 1)*(donorPts_[5] - donorPts_[4] + 1); 

		E_Int ipt_ind_dm_omp_thread[6]; 

		indice_boucle_lu_(NoD, ithread_loc, Nbre_thread_actif_loc, param_int[NoD][ ITYPCP ],
				  donorPts_,
				  topology, ipt_ind_dm_omp_thread);
		E_Int donorPts[6];
		donorPts[0]=ipt_ind_dm_omp_thread[0];
		donorPts[1]=ipt_ind_dm_omp_thread[1];
		donorPts[2]=ipt_ind_dm_omp_thread[2];
		donorPts[3]=ipt_ind_dm_omp_thread[3];
		donorPts[4]=ipt_ind_dm_omp_thread[4];
		donorPts[5]=ipt_ind_dm_omp_thread[5];
