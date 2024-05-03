# include "Fast/fast.h"
# include "Fast/param_solver.h"
# include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace K_FLD;
using namespace std;



//=============================================================================
// calcul Sij sur les points donneur des zones NS pour transfert vers zone LBM
//=============================================================================
void K_FAST::compute_sij( 
    E_Float**& ipt_ro      , E_Float**& iptS   , E_Float**& ipt_vol   , E_Int*& param_int_tc,
    E_Float*& param_real_tc, E_Int**& param_int, E_Float**& param_real, E_Int*& ipt_omp,
    E_Int& TypeTransfert   , E_Int& it_target  , E_Int& nidom         , E_Int& NoTransfert,
    E_Int& bidim           , E_Int& nstep      , E_Int& nssiter       , E_Int& rk, E_Int& exploc, E_Int& num_passage)
{

#ifdef TimeShow
  E_Float time_in = omp_get_wtime();
#endif


  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL

  // nvars pour dimensionner vectOfDnrFields
  E_Int nvars=11;

  E_Int sizecomIBC = param_int_tc[2];
  E_Int sizecomID  = param_int_tc[3+sizecomIBC];
  E_Int shift_graph = sizecomIBC + sizecomID + 3;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech           = param_int_tc[NoTransfert + shift_graph];
  E_Int nrac          = param_int_tc[ech + 1];  // nb total de raccord
  E_Int nrac_inst     = param_int_tc[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel     = param_int_tc[ech + 3];  // nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady   = nrac - nrac_inst;        // nb total de raccord stationnaire

  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;

  if (nrac_inst > 0) pass_inst_fin=2;

  //on optimise les transfert pour implicit local
  E_Int impli_local[nidom];
  if (ipt_omp == NULL) // Merde fastP a gerer
    { for (E_Int nd = 0; nd < nidom; nd++) {impli_local[nd]=1;}  }
  else{ 
        E_Int nbtask = ipt_omp[nstep-1]; 
        E_Int ptiter = ipt_omp[nssiter+ nstep-1];

        for (E_Int nd = 0; nd < nidom; nd++) {impli_local[nd]=1;}

        for (E_Int ntask = 0; ntask < nbtask; ntask++)
        {
          E_Int pttask = ptiter + ntask*(6+threadmax_sdm*7);
          E_Int nd = ipt_omp[ pttask ];
          impli_local[nd]=1;
        }
      }

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(  size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation]; // Pour l explicite local

  // on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx = 0; E_Int ibcTypeMax=0;
  for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
  {
    E_Int irac_deb = 0;
    E_Int irac_fin = nrac_steady;
    if(pass_inst == 1){ irac_deb = param_int_tc[ ech + 4 + it_target ]; irac_fin = param_int_tc[ ech + 4 + it_target + timelevel ]; }

   for (E_Int irac = irac_deb; irac < irac_fin; irac++)
    {
      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;

      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      //  si incompatibilite pass et typeTransfert, on skippe le raccord
      E_Int ibcType = param_int_tc[shift_rac + nrac * 3];

      if (ibcType > ibcTypeMax){ ibcTypeMax= ibcType;}
      E_Int ibc = 1;
      if (ibcType < 0) ibc = 0;
      if      (TypeTransfert == 0 && ibc == 1) { continue; } 
      else if (TypeTransfert == 1 && ibc == 0) { continue; }


      // Si on est en explicit local, on va autoriser les transferts entre certaines zones en fonction de la ss-ite courante
      if(exploc == 1)
	{
	  E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
	  E_Int levelD = param_int_tc[debut_rac + 25];
	  E_Int levelR = param_int_tc[debut_rac + 24];
	  E_Int cyclD  = nssiter/levelD;

	  // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	  if (levelD > levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==cyclD-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1)) { autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	  // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	  else if (levelD < levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)
		{ autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	  // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	  else if (levelD == levelR and num_passage == 1)
	    {
	      if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1)
		{ autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	  // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
	  else if (levelD == levelR and num_passage == 2)
	    {
	      if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1) { autorisation_transferts[pass_inst][irac_auto]=1; }
	    }
	}
      // Sinon, on autorise les transferts  si la zone donneuse a ete modifiee a l'iteration nstep
      else {
             E_Int NoD      =  param_int_tc[ shift_rac + nrac*5     ];
             if (impli_local[NoD]==1) autorisation_transferts[pass_inst][irac_auto]=1;
           }

      if (autorisation_transferts[pass_inst][irac_auto]==1)
	   { 
	      E_Int nbRcvPts = param_int_tc[shift_rac + nrac * 10 + 1];
    
	      if (nbRcvPts > nbRcvPts_mx) nbRcvPts_mx = nbRcvPts;

           }// autorisation transfert

    }
  }

  E_Int size = (nbRcvPts_mx / threadmax_sdm) + 1;  // on prend du gras pour gerer le residus
  E_Int r = size % 8;
  if (r != 0) size = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <= 1) size = 0;        // tableau inutile

  FldArrayF tmp(size * 17 * threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

//# pragma omp parallel default(shared)  num_threads(1)
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int ithread           = omp_get_thread_num() + 1;
    E_Int Nbre_thread_actif = omp_get_num_threads();  // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc /*, nocf*/, indCoef, noi, sizecoefs,
        /*Nbchunk,*/ imd, jmd, imdjmd;

    vector<E_Float*> vectOfDnrFields(nvars);

    // 1ere pass_typ: IBC
    // 2eme pass_typ: transfert
    //

    E_Int count_racIBC = 0;

     for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++)
     {
      // 1ere pass_inst: les raccord fixe
      // 2eme pass_inst: les raccord instationnaire
      for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
      {
        //  printf("pass %d %d %d %d \n", rank ,ipass_typ, pass_inst, it_target ); 
        // printf("ipass_inst = %d, level= %d \n",  ipass_inst, nrac_inst_level
        // );
        E_Int irac_deb = 0;
        E_Int irac_fin = nrac_steady;
        if (pass_inst == 1) {
          irac_deb = param_int_tc[ech + 4 + it_target];
          irac_fin = param_int_tc[ech + 4 + it_target + timelevel];
        }

         //printf("iracdeb=  %d, iracfin= %d \n", irac_deb, irac_fin  );
        for (E_Int irac = irac_deb; irac < irac_fin; irac++) {

	 E_Int irac_auto= irac-irac_deb;
	 if (autorisation_transferts[pass_inst][irac_auto]==1)
	  {

          E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
          //printf("ipass_typ = %d, pass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,pass_inst,irac , ithread );
          // ipass_typ,ipass_inst,irac , ithread );
          E_Int ibcType = param_int_tc[shift_rac + nrac * 3];
          E_Int ibc = 1;
          if (ibcType < 0) ibc = 0;

          E_Int NoD       = param_int_tc[shift_rac + nrac * 5     ];
          E_Int NoR       = param_int_tc[shift_rac + nrac * 11 + 1];
          E_Int nvars_loc = param_int_tc[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES

          // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
          E_Int solver_D=2; E_Int solver_R=2;
          if (nvars_loc == 11) {solver_R =4;}
          if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
          if (nvars_loc == 19) {solver_D =4; solver_R=4;}

          //printf("nvar loc %d , solver_RD= %d %d irc = %d \n", nvars_loc, solver_R, solver_D, irac);

          E_Int meshtype = param_int[ NoD ][ MESHTYPE ] ;

          if (nvars_loc == 11 ) // //Transfert NS -> LBM    
          {
            // On commence par copier les 5 variables macros
            for (E_Int eq = 0; eq < 5; eq++)
            {
              vectOfDnrFields[eq] = ipt_ro[NoD] + eq * param_int[NoD][ NDIMDX];
            }
            // Puis on copie les gradients
            for (E_Int eq = 5; eq < nvars_loc; eq++)
            {
              vectOfDnrFields[eq] = iptS[NoD] + (eq-5) * param_int[NoD][ NDIMDX];
            }
          }
          else { continue;}


          imd = param_int[ NoD ][NIJK  ];
          jmd = param_int[ NoD ][NIJK+1];
          // }

          imdjmd = imd * jmd;

          ////
          //  Interpolation parallele
          ////
          ////

          E_Int nbRcvPts = param_int_tc[shift_rac + nrac * 10 + 1];
          // E_Int nbDonPts = param_int_tc[ shift_rac                ];

          E_Int pos;
          pos = param_int_tc[shift_rac + nrac * 7];      E_Int* ntype    = param_int_tc  + pos;
          pos = pos + 1 + ntype[0];                      E_Int* types    = param_int_tc  + pos;
          pos = param_int_tc[shift_rac + nrac * 6];      E_Int* donorPts = param_int_tc  + pos;
          pos = param_int_tc[shift_rac + nrac * 12 + 1]; E_Int* rcvPts   = param_int_tc  + pos;  // donor et receveur inverser car storage donor
          pos = param_int_tc[shift_rac + nrac * 8];    E_Float* ptrCoefs = param_real_tc + pos;

          E_Int nbInterpD = param_int_tc[shift_rac + nrac];

          E_Int ideb = 0;
          E_Int ifin = 0;
          E_Int shiftCoef = 0;
          E_Int shiftDonor = 0;

          for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++) {
            type = types[ifin];

            SIZECF(type, meshtype, sizecoefs);
            ifin = ifin + ntype[1 + ndtyp];


            E_Int pt_deb, pt_fin;

            /// oldschool
            // Calcul du nombre de champs a traiter par chaque thread
            E_Int size_bc = ifin - ideb;
            E_Int chunk = size_bc / Nbre_thread_actif;
            E_Int r = size_bc - chunk * Nbre_thread_actif;
            // pts traitees par thread
            if (ithread <= r) {
              pt_deb = ideb + (ithread - 1) * (chunk + 1);
              pt_fin = pt_deb + (chunk + 1);
            } else {
              pt_deb = ideb + (chunk + 1) * r + (ithread - r - 1) * chunk;
              pt_fin = pt_deb + chunk;
            }

            // Si type 0, calcul sequentiel
            if (type == 0) {
              if (ithread == 1) {
                pt_deb = ideb;
                pt_fin = ifin;
              } else {
                pt_deb = ideb;
                pt_fin = ideb;
              }
            }

                      noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
                      indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
                      E_Int shiftv =0;

                      if (solver_D<4 && solver_R==4)
                      {
                         // Transfert NS vers LBM : adimensionnement
#                        include "Fast/INTERP/cp_sij.h"
                      }


                      ideb       = ideb + ifin;
                      shiftCoef  = shiftCoef +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
                      shiftDonor= shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif

                   }// type
	          } //autorisation transfert
                }//irac
               }//pass_inst
    }  // ipass
  }    // omp

}
