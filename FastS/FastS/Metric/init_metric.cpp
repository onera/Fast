/*    
    Copyright 2013-2022 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
 
#include "fastS.h"
#include "param_solver.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
using namespace K_FLD;
#include "stub.h"
extern int __activation__;

//=============================================================================
// La metrique
// IN: arrays: les coords des zones
// IN: num: le dictionnaires numerics
//=============================================================================
PyObject* K_FASTS::init_metric(PyObject* self, PyObject* args)
{
  if (__activation__ == 0) { PyErr_SetString(PyExc_NotImplementedError, STUBMSG); return NULL; }

  PyObject* zones;  PyObject* metrics; E_Int omp_mode;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOl", &zones, &metrics, &omp_mode)) return NULL; 
#else
  if (!PyArg_ParseTuple(args, "OOi", &zones, &metrics, &omp_mode)) return NULL; 
#endif

  E_Int nidom    = PyList_Size(zones);

  E_Int** ipt_ind_dm; E_Int** ipt_degen;
  
  E_Float** param_real;
  E_Float** ipty;       E_Float** iptz;  
  E_Float** ipti;       E_Float** iptj;     E_Float** iptk;    E_Float** iptvol;
  E_Float** ipti_df;    E_Float** iptj_df;  E_Float** iptk_df; E_Float** iptvol_df;
  E_Float** ipti0;      E_Float** iptj0;    E_Float** iptk0;   
  E_Float** iptventi;   E_Float** iptventj; E_Float** iptventk; E_Float** iptdist;

  E_Int* param_int[nidom*3];

  ipt_ind_dm        = param_int   + nidom;
  ipt_degen         = ipt_ind_dm      + nidom;

  E_Float* iptx[nidom*19];
  ipty              = iptx            + nidom;
  iptz              = ipty            + nidom;
  ipti              = iptz            + nidom;
  iptj              = ipti            + nidom;
  iptk              = iptj            + nidom;
  iptvol            = iptk            + nidom;
  ipti0             = iptvol          + nidom;
  iptj0             = ipti0           + nidom;
  iptk0             = iptj0           + nidom;
  ipti_df           = iptk0           + nidom;
  iptj_df           = ipti_df         + nidom;
  iptk_df           = iptj_df         + nidom;
  iptvol_df         = iptk_df         + nidom;
  iptventi          = iptvol_df       + nidom;
  iptventj          = iptventi        + nidom;
  iptventk          = iptventj        + nidom;
  param_real    = iptventk        + nidom;
  iptdist           = param_real  + nidom;

  vector<PyArrayObject*> hook;



  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    // check zone
    PyObject* zone = PyList_GetItem(zones, nd); // domaine i


    /* Get numerics from zone */
    PyObject* numerics    = K_PYTREE::getNodeFromName1(zone, ".Solver#ownData");
    PyObject*          t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    param_int[nd]     = K_PYTREE::getValueAI(t, hook);
                       t  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    param_real[nd]    = K_PYTREE::getValueAF(t, hook);

    /*-------------------------------------*/
    /* Extraction (x,y,z): pour forcage spatia */
    /*-------------------------------------*/
    if(param_int[nd][ LALE ]== 0 || param_int[nd][ LALE ]== 2 ){ GET_XYZ( "GridCoordinates"     , zone, iptx[nd], ipty[nd], iptz[nd])}
    else                                                       { GET_XYZ( "GridCoordinates#Init", zone, iptx[nd], ipty[nd], iptz[nd])}
    //GET_XYZ( "GridCoordinates", zone, iptx[nd], ipty[nd], iptz[nd]) 

    /* get distance paroi */

    //Pointeur visqueux: mut, dist, zgris sont en acces compact
    if(param_int[nd][ IFLOW ] ==  3)
      { PyObject* sol_center;
        sol_center  = K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
        t           = K_PYTREE::getNodeFromName1(sol_center, "TurbulentDistance");
        iptdist[nd] = K_PYTREE::getValueAF(t, hook);
      }
    else {iptdist[nd] = NULL;}

    /* get metric */
    PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i

    E_Float* dummy;
    E_Int lale, kfludom;
    lale = param_int[nd][LALE]; kfludom = param_int[nd][KFLUDOM];
    GET_TI(METRIC_TI  , metric, param_int[nd], ipti [nd]  , iptj [nd]  , iptk [nd]  , iptvol[nd]   )
    if (lale>=1)
     {
       GET_TI(METRIC_TI0 , metric, param_int[nd], ipti0[nd]  , iptj0[nd]  , iptk0[nd]  , dummy )
       GET_VENT( METRIC_VENT, metric, param_int[nd], iptventi[nd], iptventj[nd], iptventk[nd] )
     }
    else
     {   ipti0[nd]    = ipti[nd]; iptj0[nd]=iptj[nd]; iptk0[nd]=iptk[nd]; 
         iptventi[nd] = ipti[nd]; iptventj[nd] = ipti[nd]; iptventk[nd] = ipti[nd];
     }
    if(kfludom == 3)  
     {
       GET_TI(METRIC_TIDF, metric, param_int[nd], ipti_df[nd], iptj_df[nd], iptk_df[nd], iptvol_df[nd])
     }
    else
     {   ipti_df[nd]    = ipti[nd]; iptj_df[nd]=iptj[nd]; iptk_df[nd]=iptk[nd]; 
     }

    ipt_degen[nd ]   = K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_DEGEN) );
    ipt_ind_dm[ nd ] = K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM)  );

  } // boucle zone


  //* tableau pour stocker dimension sous-domaine omp *//
  E_Int threadmax_sdm = 1;
#ifdef _OPENMP
  threadmax_sdm = omp_get_max_threads();
#endif
  
  FldArrayI ind_sdm(6*threadmax_sdm);          E_Int* ipt_ind_sdm         = ind_sdm.begin();
  FldArrayI ind_coe(6*threadmax_sdm);          E_Int* ipt_ind_coe         = ind_coe.begin();
  FldArrayI ind_grad(6*threadmax_sdm);         E_Int* ipt_ind_grad        = ind_grad.begin();
  FldArrayI ind_dm1(6*threadmax_sdm);          E_Int* ipt_ind_dm1         = ind_dm1.begin();
  FldArrayI ijkv_sdm(3*threadmax_sdm);         E_Int* ipt_ijkv_sdm        = ijkv_sdm.begin();
  FldArrayI ind_dm_omp(12*threadmax_sdm);      E_Int* ipt_ind_dm_omp      = ind_dm_omp.begin();
  FldArrayI topology_socket(3*threadmax_sdm);  E_Int* ipt_topology_socket = topology_socket.begin();

  //FldArrayF rot_ale(12*threadmax_sdm);         E_Float* ipt_rot_ale       = rot_ale.begin();

#pragma omp parallel default(shared) 
     {
	//* variable declaree dans zone parallele = private *//
#ifdef _OPENMP
        E_Int  ithread           = omp_get_thread_num() +1;
        E_Int  Nbre_thread_actif = omp_get_num_threads();
#else
        E_Int  ithread           = 1;
        E_Int  Nbre_thread_actif = 1;
#endif
        //E_Int Nbre_socket   = NBR_SOCKET;                       // nombre de proc (socket) sur le noeud a memoire partagee
        E_Int Nbre_socket   = 1;                       // nombre de proc (socket) sur le noeud a memoire partagee
        if( Nbre_thread_actif < Nbre_socket ) Nbre_socket = 1;

        E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
        E_Int socket          = (ithread-1)/thread_parsock +1;
        E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

        E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm    + (ithread-1)*3; 
        E_Int* ipt_ind_sdm_thread     = ipt_ind_sdm     + (ithread-1)*6;
        E_Int* ipt_ind_coe_thread     = ipt_ind_coe     + (ithread-1)*6;
        E_Int* ipt_ind_grad_thread    = ipt_ind_grad    + (ithread-1)*6;

        E_Int* ipt_ind_dm_loc        = ipt_ind_dm1     + (ithread-1)*6;
        E_Int* ipt_ind_dm_socket     = ipt_ind_dm_omp  + (ithread-1)*12;
        E_Int* ipt_ind_dm_omp_thread = ipt_ind_dm_socket  + 6;

        E_Int* ipt_topology_socket_thread =  ipt_topology_socket + (ithread-1)*3;

        E_Int* ipt_topo_omp; E_Int* ipt_inddm_omp; E_Int ithread_loc; E_Int Nbre_thread_actif_loc;
        if (omp_mode == 1)
            { 
              // loop calcul normale
              for (E_Int nd = 0; nd < nidom; nd++)
              {
#              include "Metric/indice_omp1.h" 
               cp_tijk_( param_int[nd], iptx[nd], ipty[nd], iptz[nd], ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], ind_mtr);
              }
       	      #pragma omp barrier
              // loop calcul volume
              for (E_Int nd = 0; nd < nidom; nd++)
              {
#              include "Metric/indice_omp1.h" 
               cp_vol_(  param_int[nd], iptx[nd], ipty[nd], iptz[nd], ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], iptvol[nd], ind_mtr);

               for (E_Int k = ind_mtr[4]; k <= ind_mtr[5]; k++){ 
                for (E_Int j = ind_mtr[2]; j <= ind_mtr[3]; j++){ 
                 for (E_Int i = ind_mtr[0]; i <= ind_mtr[1]; i++){ 

                   E_Int l =  (i+ param_int[nd][NIJK_MTR+3]-1)*param_int[nd][NIJK_MTR]
                            + (j+ param_int[nd][NIJK_MTR+3]-1)*param_int[nd][NIJK_MTR+1]
                            + (k+ param_int[nd][NIJK_MTR+4]-1)*param_int[nd][NIJK_MTR+2];
                   iptvol[nd][l] = K_FUNC::E_max(iptvol[nd][l], 1.e-30);
                 }
                }
               }
              }
       	      #pragma omp barrier
              // loop extrapolation
              for (E_Int nd = 0; nd < nidom; nd++)
              {

#              include "Metric/indice_omp1.h" 
               if (param_int[nd][ ITYPZONE ] != 2  && ithread_loc == 1)
                {
                 ipt_ind_dm_loc[0]   = 1; 
                 ipt_ind_dm_loc[2]   = 1; 
                 ipt_ind_dm_loc[4]   = 1; 
                 if(param_int[nd][ ITYPZONE ] == 0) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = param_int[nd][ IJKV+2 ];
                   }
                 else if(param_int[nd][ ITYPZONE ] == 1) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1;
                   }
                 else if(param_int[nd][ ITYPZONE ] == 2) 
                   {
                    ipt_ind_dm_loc[1]   = 1 ;
                    ipt_ind_dm_loc[3]   = 1; 
                    ipt_ind_dm_loc[5]   = 1; 
                   }
                 else
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1 ;
                   }

                 E_Int* ipt_nijk_xyz = param_int[nd]+ NIJK_XYZ;
                 E_Int* ipt_nijk_mtr = param_int[nd]+ NIJK_MTR;
                 E_Int* ipt_nijk     = param_int[nd]+ NIJK;

                 tijk_extrap_( param_int[nd][ NDIMDX_MTR ], param_int[nd][ NDIMDX_XYZ ] , ipt_nijk_xyz, ipt_nijk_mtr,
                               param_int[nd][ NEQ_IJ ]    , param_int[nd][ NEQ_K ],
                               ipt_ind_dm_loc,
                               ipt_degen[nd] ,
                               ipti[nd]      , iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], iptvol[nd]); 

                 if(param_int[nd][ IFLOW ] ==3)
                   { ipt_ind_dm_loc[1]= param_int[nd][ IJKV ]; ipt_ind_dm_loc[3]= param_int[nd][ IJKV+1 ]; ipt_ind_dm_loc[5]= param_int[nd][IJKV+2];

                     dist_extrap_( param_int[nd][ NDIMDX ], param_int[nd][ NDIMDX_XYZ ] , ipt_nijk, ipt_nijk_xyz,
                                    ipt_ind_dm_loc, ipt_degen[nd] , iptdist[nd]);
                   }
                }
              }
            }
        else
            { //omp_mode 0
              for (E_Int nd = 0; nd < nidom; nd++)
              {
                E_Int lmin = 10;
                if (param_int[nd][ ITYPCP ] == 2) lmin = 4;

                ipt_ind_dm_loc[0]   = 1; 
                ipt_ind_dm_loc[2]   = 1; 
                ipt_ind_dm_loc[4]   = 1; 
                if(param_int[nd][ ITYPZONE ] == 0) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = param_int[nd][ IJKV+2 ];
                   }
                else if(param_int[nd][ ITYPZONE ] == 1) 
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1;
                   }
                else if(param_int[nd][ ITYPZONE ] == 2) 
                   {
                    ipt_ind_dm_loc[1]   = 1 ;
                    ipt_ind_dm_loc[3]   = 1; 
                    ipt_ind_dm_loc[5]   = 1; 
                   }
                else
                   {
                    ipt_ind_dm_loc[1]   = param_int[nd][ IJKV   ];
                    ipt_ind_dm_loc[3]   = param_int[nd][ IJKV+1 ];
                    ipt_ind_dm_loc[5]   = 1 ;
                   }

                 indice_boucle_lu_(nd, socket , Nbre_socket, lmin,
                                   ipt_ind_dm_loc, 
                                   ipt_topology_socket_thread, ipt_ind_dm_socket );

                 skmtr_( nd, param_int[nd], param_real[nd],
	                iptx[nd], ipty[nd], iptz[nd], ipt_degen[nd], iptdist[nd],
                        ipti[nd], iptj[nd], iptk[nd], ipti0[nd], iptj0[nd], iptk0[nd], iptvol[nd], iptventi[nd], iptventj[nd], iptventk[nd],
                        ipt_ijkv_sdm_thread,
                        ipt_ind_sdm_thread , ipt_ind_coe_thread, ipt_ind_grad_thread    , 
                        ipt_ind_dm_loc     , ipt_ind_dm_socket , ipt_ind_dm_omp_thread  ,
                        ipt_topology_socket_thread, 
                        ithread_sock       , thread_parsock    , Nbre_thread_actif, Nbre_socket, socket,
                        ithread);

                 // cutoff vol
                 #pragma omp for
                 for (E_Int i = 0; i < param_int[nd][NDIMDX_MTR]; i++) { iptvol[nd][i] = K_FUNC::E_max(iptvol[nd][i], 1.e-30);}

	         if( param_int[nd][KFLUDOM] == 3)
	         {

       	          #pragma omp single
     	          {
                  //* tableau pour calculer dimension metric DF du domaine *//
                  FldArrayF  iptmpi(   param_int[nd][ NDIMDX_XYZ ]);
                  FldArrayF  iptmpj(   param_int[nd][ NDIMDX_XYZ ]);
                  FldArrayF  iptmpk(   param_int[nd][ NDIMDX_XYZ ]);
                  FldArrayF iptmpi2(   param_int[nd][ NDIMDX_XYZ ]);
                  FldArrayF iptmpj2(   param_int[nd][ NDIMDX_XYZ ]);
                  FldArrayF iptmpk2(   param_int[nd][ NDIMDX_XYZ ]);
                  FldArrayF  iptmtr( 9*param_int[nd][ NDIMDX_XYZ ]);

                  E_Float* tmpi = iptmpi.begin();
                  E_Float* tmpj = iptmpj.begin();
                  E_Float* tmpk = iptmpk.begin();
                  E_Float* tmpi2= iptmpi2.begin();
                  E_Float* tmpj2= iptmpj2.begin();
                  E_Float* tmpk2= iptmpk2.begin();
                  E_Float* mtr  = iptmtr.begin();

                  skmtr_df_( param_int[nd]+NIJK_XYZ, param_int[nd]+NIJK_MTR,
	                     param_int[nd][ NDIMDX_XYZ ], iptx[nd], ipty[nd], iptz[nd], 
	                     tmpi,tmpj,tmpk,tmpi2,tmpj2,tmpk2, mtr,
	                     param_int[nd][ NEQ_IJ ], param_int[nd][ NEQ_K ],
		             param_int[nd][ NDIMDX_MTR], ipti_df[nd],iptj_df[nd],iptk_df[nd], iptvol_df[nd]);
                  } //omp single
	         } // if kflu=DF
              }//zone
           }// omp_mode


     } //* fin zone parallele  *//

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
 
  Py_INCREF(Py_None);
  return Py_None;
}
