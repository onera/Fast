/*    
    Copyright 2013-2019 Onera.

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
PyObject* K_FASTS::allocate_metric(PyObject* self, PyObject* args)
{
  if (__activation__ == 0) { PyErr_SetString(PyExc_NotImplementedError, STUBMSG); return NULL; }

  PyObject *zone; E_Int nssiter;
#if defined E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &zone, &nssiter)) return NULL; 
#else 
  if (!PyArg_ParseTuple(args, "Oi", &zone, &nssiter)) return NULL; 
#endif

  vector<PyArrayObject*> hook;
  PyObject* metrics  = PyList_New(0); 

  E_Float* iptx; E_Float* ipty; E_Float* iptz;

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

  FldArrayF rot_ale(12*threadmax_sdm);         E_Float* ipt_rot_ale       = rot_ale.begin();

  //printf("Nombre de thread max= %d\n", threadmax_sdm);


    /* Get numerics from zone */
    PyObject*   numerics    = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*            o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    E_Int* ipt_param_int    = K_PYTREE::getValueAI(o, hook);
                         o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    E_Float* ipt_param_real = K_PYTREE::getValueAF(o, hook);

    E_Int lale, kfludom;
    lale = ipt_param_int[LALE]; kfludom = ipt_param_int[KFLUDOM]; E_Int* size_ssdom = ipt_param_int+ SIZE_SSDOM;

    //
    //
    //Pointeur maillage
    //
    //
    PyObject* sol_center; PyObject* t; PyObject* ro;
    if(lale==1){sol_center   = K_PYTREE::getNodeFromName1(zone      , "GridCoordinates#Init");}
    else       {sol_center   = K_PYTREE::getNodeFromName1(zone      , "GridCoordinates"     );}

    if (sol_center == NULL) { PyErr_SetString(PyExc_ValueError, "Coord missing for metric computation"); return NULL; }

    t        = K_PYTREE::getNodeFromName1(sol_center, "CoordinateX");
    iptx     = K_PYTREE::getValueAF(t, hook);
    t        = K_PYTREE::getNodeFromName1(sol_center, "CoordinateY");
    ipty     = K_PYTREE::getValueAF(t, hook);
    t        = K_PYTREE::getNodeFromName1(sol_center, "CoordinateZ");
    iptz     = K_PYTREE::getValueAF(t, hook);

    sol_center= K_PYTREE::getNodeFromName1(zone      , "FlowSolution#Centers");
    if (sol_center == NULL) ro = t;                                        //blindage pour appeler metric sur un arbre qui contient seuleemt x,y,x
    else
      {t         = K_PYTREE::getNodeFromName1(sol_center, "Density");      //Objet python= noeud CGNS
       ro        = PyList_GetItem(t, 1);                                   //Objet python= tableau numpy du noeud CGNS
      }

    //
    // On recupere taille zone et dimensionne tableau
    //
    E_Int* dim = K_PYTREE::getValueAI(zone, hook);
    E_Int ni =  dim[0]; E_Int nj = dim[1]; E_Int nk = dim[2];

    E_Int nfic =2;
    E_Int ific = nfic; E_Int ific_xyz = nfic;
    if (kfludom == 3)  { ific_xyz = 6; ific = 3; }

    //* tableau pour stocker dimension x,y,z du domaine *//
    ipt_param_int[ NIJK_XYZ   ]  = ni;
    ipt_param_int[ NIJK_XYZ +1]  = nj;
    ipt_param_int[ NIJK_XYZ +2]  = nk;
    ipt_param_int[ NIJK_XYZ +3]  = ific_xyz;
    ipt_param_int[ NIJK_XYZ +4]  = ific_xyz;
    if (nk == 2) ipt_param_int[ NIJK_XYZ +4]=0;

    //* tableau pour stocker rop du domaine *//
    ipt_param_int[ NIJK   ]  = ni-1;
    ipt_param_int[ NIJK +1]  = nj-1;
    ipt_param_int[ NIJK +2]  = nk-1;
    ipt_param_int[ NIJK +3]  = ific;
    ipt_param_int[ NIJK +4]  = ific;
    if (nk == 2) ipt_param_int[ NIJK +4]= 0;

    //* tableau pour stocker dimension rop du domaine sans Ghostcell*//
    ipt_param_int[ IJKV   ]  = ipt_param_int[ NIJK   ]-2*ipt_param_int[ NIJK +3];
    ipt_param_int[ IJKV +1]  = ipt_param_int[ NIJK+1 ]-2*ipt_param_int[ NIJK +3];
    ipt_param_int[ IJKV +2]  = ipt_param_int[ NIJK+2 ]-2*ipt_param_int[ NIJK +4];
    
    ipt_param_int[ NDIMDX      ] =  ipt_param_int[ NIJK     ]*ipt_param_int[ NIJK+1     ]*ipt_param_int[ NIJK+2     ] + ipt_param_int[SHIFTVAR]; 
    ipt_param_int[ NDIMDX_XYZ  ] =  ipt_param_int[ NIJK_XYZ ]*ipt_param_int[ NIJK_XYZ+1 ]*ipt_param_int[ NIJK_XYZ+2 ];

    //
    //* determine la nature de la zone: 3d curvi, 3d-kcart, 2d,...  *//
    //
    E_Int neq_ij = 0; E_Int neq_k = 0; 
    E_Int typ_zone, ndimdx_mtr ; 

    
    //FldArrayI  degener(ipt_param_int[ NDIMDX_XYZ ]); E_Int* ipt_degener = degener.begin();
    PyObject* degener    = K_NUMPY::buildNumpyArray(ipt_param_int[ NDIMDX_XYZ ], 1, 1, 1);
    E_Int*   ipt_degener = K_NUMPY::getNumpyPtrI( degener );

    nature_geom_dom_(ipt_param_int+NIJK_XYZ, ipt_param_int[ NDIMDX_XYZ ], iptx, ipty, iptz , ipt_degener, lale, ipt_param_real, // IN
                     ipt_param_int+NIJK_MTR, ndimdx_mtr, neq_ij, neq_k, typ_zone);                             // OUT

    ipt_param_int[ NDIMDX_MTR  ]  = ndimdx_mtr;
    ipt_param_int[ ITYPZONE    ]  = typ_zone;
    ipt_param_int[ NEQ_IJ      ]  = neq_ij;
    ipt_param_int[ NEQ_K       ]  = neq_k;

    E_Int coe_shift = 0;
    //determination du nombre de "sous-pas/sousiter" en fonction du domaine impli/expli
    if     (ipt_param_int[ ITYPCP ] <= 1) { coe_shift = 0; }
    else if(ipt_param_int[ ITYPCP ] == 2) { coe_shift = 1; }
    else { printf("Warning: invalid temporal scheme."); }
          
    if  (ipt_param_int[ IFLOW ] == 3 && ipt_param_int[ ILES ] == 0){  ipt_param_int[ NEQ_COE ] = 6 -coe_shift*5;}    //!6 en implicite, 1 en rk3
    else{                                                             ipt_param_int[ NEQ_COE ] = 5 -coe_shift*4;}    //!5 en implicite, 1 en rk3
    
    if  (ipt_param_int[ IFLOW ] == 3){ ipt_param_int[ NEQ ] = 6;} 
    else                             { ipt_param_int[ NEQ ] = 5;}

    //determination taille vitesse entrainement
     E_Float rot_z, rot_tot;
    if(lale==1)
    {
      rot_z  =  ipt_param_real[ROT_OMEGA  ]*ipt_param_real[ROT_OMEGA ]+ ipt_param_real[ROT_OMEGA+1 ]*ipt_param_real[ROT_OMEGA+1 ];
      rot_tot=  rot_z + ipt_param_real[ROT_OMEGA+2 ]*ipt_param_real[ROT_OMEGA+2 ];
    }
    else
    { rot_z  = 0.; rot_tot= 0.; }

    ipt_param_int[ ITYPVENT ] = 0;
    if  (rot_z   < 1.e-14)             ipt_param_int[ ITYPVENT ] = 1; //rot axe z
    if  (rot_tot < 1.e-14)             ipt_param_int[ ITYPVENT ] = 2; //translation
    if  (ipt_param_int[ ITYPZONE ]==3) ipt_param_int[ ITYPVENT ] = 3; // maillage 2d

     //printf("vent type= %d  \n", ipt_param_int[ ITYPVENT ]);

    if      (ipt_param_int[ ITYPVENT ] == 0)
        { 
           ipt_param_int[ NEQ_VENT ]    = 3;

           ipt_param_int[ NIJK_VENT  ]  = 1;
           ipt_param_int[ NIJK_VENT+1]  = ipt_param_int[ NIJK_XYZ  ];
           ipt_param_int[ NIJK_VENT+2]  = ipt_param_int[ NIJK_XYZ  ]*ipt_param_int[ NIJK_XYZ+1 ];
           ipt_param_int[ NIJK_VENT+3]  = ipt_param_int[ NIJK_XYZ+3];
           ipt_param_int[ NIJK_VENT+4]  = ipt_param_int[ NIJK_XYZ+4];
           ipt_param_int[ NDIMDX_VENT]  = ipt_param_int[ NDIMDX_XYZ ];
        }
    else if (ipt_param_int[ ITYPVENT ] == 1)
       { 
           ipt_param_int[ NEQ_VENT ]    = 2;

           if      (ipt_param_int[ ITYPZONE ] == 0)
                {
                 ipt_param_int[ NDIMDX_VENT]  = ipt_param_int[ NDIMDX_XYZ ];

                 ipt_param_int[ NIJK_VENT  ]  = 1;
                 ipt_param_int[ NIJK_VENT+1]  = ipt_param_int[ NIJK_XYZ  ];
                 ipt_param_int[ NIJK_VENT+2]  = ipt_param_int[ NIJK_XYZ  ]*ipt_param_int[ NIJK_XYZ+1 ];
                 ipt_param_int[ NIJK_VENT+3]  = ipt_param_int[ NIJK_XYZ+3];
                 ipt_param_int[ NIJK_VENT+4]  = ipt_param_int[ NIJK_XYZ+4];
                }
           else if (ipt_param_int[ ITYPZONE ] == 1)
                {
                 ipt_param_int[ NDIMDX_VENT]  = ipt_param_int[ NIJK_XYZ  ]*ipt_param_int[ NIJK_XYZ+1 ];
                 ipt_param_int[ NIJK_VENT  ]  = 1;
                 ipt_param_int[ NIJK_VENT+1]  = ipt_param_int[ NIJK_XYZ  ];
                 ipt_param_int[ NIJK_VENT+2]  = 0 ;
                 ipt_param_int[ NIJK_VENT+3]  = ipt_param_int[ NIJK_XYZ+3];
                 ipt_param_int[ NIJK_VENT+4]  = 0;
                }
       } 
    else if (ipt_param_int[ ITYPVENT ] == 2)
       {
           ipt_param_int[ NEQ_VENT   ]  = 3;

           ipt_param_int[ NDIMDX_VENT]  = 1 ;
           ipt_param_int[ NIJK_VENT  ]  = 0;
           ipt_param_int[ NIJK_VENT+1]  = 0;
           ipt_param_int[ NIJK_VENT+2]  = 0 ;
           ipt_param_int[ NIJK_VENT+3]  = 0;
           ipt_param_int[ NIJK_VENT+4]  = 0;

       }
    else 
       {
           ipt_param_int[ NEQ_VENT ]    = 2;

           ipt_param_int[ NDIMDX_VENT]  = ipt_param_int[ NIJK_XYZ  ]*ipt_param_int[ NIJK_XYZ+1 ];
           ipt_param_int[ NIJK_VENT  ]  = 1;
           ipt_param_int[ NIJK_VENT+1]  = ipt_param_int[ NIJK_XYZ  ];
           ipt_param_int[ NIJK_VENT+2]  = 0 ;
           ipt_param_int[ NIJK_VENT+3]  = ipt_param_int[ NIJK_XYZ+3];
           ipt_param_int[ NIJK_VENT+4]  = ipt_param_int[ NIJK_XYZ+4];
       }

     // On remonte la valeur de type zone dans l'arbre
     o = K_PYTREE::getNodeFromName1(numerics, "type_zone");
     if (o == NULL) { PyErr_SetString(PyExc_ValueError, "metric: type zone is missing or is invalid."); return 0; }
     E_Int* d = K_PYTREE::getValueAI(o, hook); d[0] = typ_zone;
     PyObject* zname_py = PyList_GetItem(zone, 0); 
     char* zname;
     if (PyString_Check(zname_py)) zname = PyString_AsString(zname_py);
#if PY_VERSION_HEX >= 0x03000000
     else if (PyUnicode_Check(zname_py)) zname = PyUnicode_AsUTF8(zname_py);
#endif
     else zname = NULL;
     if      (typ_zone == 0) {printf("typezone: 3D curvilinear, %s (%d, %d, %d)\n", zname, ipt_param_int[ IJKV ], ipt_param_int[ IJKV +1], ipt_param_int[ IJKV +2]);}
     else if (typ_zone == 1) {printf("typezone: 3D, homogenous k direction with constant step, %s (%d, %d, %d)\n", zname,ipt_param_int[ IJKV ], ipt_param_int[ IJKV +1],  ipt_param_int[ IJKV +2]);}
     else if (typ_zone == 2) {printf("typezone: 3D cartesian with constant step, %s (%d, %d, %d)\n",  zname,ipt_param_int[ IJKV ],  ipt_param_int[ IJKV +1],  ipt_param_int[ IJKV +2]);}
     else if (typ_zone == 3) {printf("typezone: 2D curvilinear, %s (%d, %d, %d)\n", zname, ipt_param_int[ IJKV ], ipt_param_int[ IJKV +1], ipt_param_int[ IJKV +2]);}

     //
     //* Declare memoire pour metric: normales + volume)
     //
     PyObject* ipti;

     E_Int neq_mtr = 2*neq_ij+ neq_k + 1; // ti+ tj+ tk+ vol
     ipti  = K_NUMPY::buildNumpyArray( ndimdx_mtr, neq_mtr, 0, 1);  
     //
     //* Declare memoire pour metric instant initial si zone ale: normales + volume)
     //
     PyObject* ipti0; PyObject* iptventi;
     E_Int neq_vent;
     if     (lale == 0) { ipti0 = ipti; iptventi= ro; neq_vent = 0; }
     else // ale=1
     {
       neq_vent = ipt_param_int[ NEQ_VENT ]*3;
       if      (ipt_param_int[ ITYPVENT ] == 2) neq_vent = ipt_param_int[ NEQ_VENT ];
       else if (ipt_param_int[ ITYPVENT ] == 3) neq_vent = ipt_param_int[ NEQ_VENT ]*2;
                     
       ipti0    = K_NUMPY::buildNumpyArray( ipt_param_int[ NDIMDX_MTR  ], neq_mtr-1 , 0, 1);  
       iptventi = K_NUMPY::buildNumpyArray( ipt_param_int[ NDIMDX_VENT ], neq_vent  , 0, 1); 
     }
     //
     //* Declare memoire pour metric DF
     //
     PyObject* ipti_df; 
     if(kfludom == 3)  ipti_df = K_NUMPY::buildNumpyArray( ndimdx_mtr, neq_mtr, 0, 1);  
     else              ipti_df = ipti; 

     E_Float* ti  = K_NUMPY::getNumpyPtrF(ipti  );
     E_Float* tj  = ti + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_IJ];
     E_Float* tk  = tj + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_IJ];
     E_Float* vol = tk + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_K ];

     E_Float* ti_df  = K_NUMPY::getNumpyPtrF(ipti_df  );
     E_Float* tj_df  = ti_df + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_IJ];
     E_Float* tk_df  = tj_df + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_IJ];
     E_Float* vol_df = tk_df + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_K ];

     E_Float* ti0  = K_NUMPY::getNumpyPtrF(ipti0  );
     E_Float* tj0  = ti0 + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_IJ];
     E_Float* tk0  = tj0 + ipt_param_int[ NDIMDX_MTR]* ipt_param_int[ NEQ_IJ];

     E_Float* venti = K_NUMPY::getNumpyPtrF(iptventi);
     E_Float* ventj = venti + ipt_param_int[ NDIMDX_VENT]*ipt_param_int[ NEQ_VENT ];
     E_Float* ventk = ventj + ipt_param_int[ NDIMDX_VENT]*ipt_param_int[ NEQ_VENT ];

     //* tableau pour calculer dimension metric DF du domaine *//
     FldArrayF  iptmpi(   ipt_param_int[ NDIMDX_XYZ ]);
     FldArrayF  iptmpj(   ipt_param_int[ NDIMDX_XYZ ]);
     FldArrayF  iptmpk(   ipt_param_int[ NDIMDX_XYZ ]);
     FldArrayF iptmpi2(   ipt_param_int[ NDIMDX_XYZ ]);
     FldArrayF iptmpj2(   ipt_param_int[ NDIMDX_XYZ ]);
     FldArrayF iptmpk2(   ipt_param_int[ NDIMDX_XYZ ]);
     FldArrayF  iptmtr( 9*ipt_param_int[ NDIMDX_XYZ ]);
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
        if (Nbre_thread_actif < Nbre_socket) Nbre_socket = 1;

        //E_Int thread_parsock  =  Nbre_thread_actif/Nbre_socket;
        //E_Int socket          = (ithread-1)/thread_parsock +1;
        //E_Int  ithread_sock   = ithread-(socket-1)*thread_parsock;

        //E_Int* ipt_ijkv_sdm_thread    = ipt_ijkv_sdm    + (ithread-1)*3; 
        //E_Int* ipt_ind_sdm_thread     = ipt_ind_sdm     + (ithread-1)*6;
        //E_Int* ipt_ind_coe_thread     = ipt_ind_coe     + (ithread-1)*6;
        //E_Int* ipt_ind_grad_thread    = ipt_ind_grad    + (ithread-1)*6;

        E_Int* ipt_ind_dm_loc        = ipt_ind_dm1     + (ithread-1)*6;
        //E_Int* ipt_ind_dm_socket     = ipt_ind_dm_omp  + (ithread-1)*12;
        //E_Int* ipt_ind_dm_omp_thread = ipt_ind_dm_socket  + 6;

        //E_Float*  ipt_rot_ale_thread = ipt_rot_ale  + (ithread-1)*12;

        ipt_ind_dm_loc[0]   = 1; 
        ipt_ind_dm_loc[2]   = 1; 
        ipt_ind_dm_loc[4]   = 1; 
        if(ipt_param_int[ ITYPZONE ] == 0) 
          {
           ipt_ind_dm_loc[1]   = ipt_param_int[ IJKV   ];
           ipt_ind_dm_loc[3]   = ipt_param_int[ IJKV+1 ];
           ipt_ind_dm_loc[5]   = ipt_param_int[ IJKV+2 ];
          }
         else if(ipt_param_int[ ITYPZONE ] == 1) 
          {
           ipt_ind_dm_loc[1]   = ipt_param_int[ IJKV   ];
           ipt_ind_dm_loc[3]   = ipt_param_int[ IJKV+1 ];
           ipt_ind_dm_loc[5]   = 1;
          }
         else if(ipt_param_int[ ITYPZONE ] == 2) 
          {
           ipt_ind_dm_loc[1]   = 1 ;
           ipt_ind_dm_loc[3]   = 1; 
           ipt_ind_dm_loc[5]   = 1; 
          }
         else
          {
           ipt_ind_dm_loc[1]   = ipt_param_int[ IJKV   ];
           ipt_ind_dm_loc[3]   = ipt_param_int[ IJKV+1 ];
           ipt_ind_dm_loc[5]   = 1 ;
          }

         E_Int nd =0;
         E_Int lmin = 10;
         if (ipt_param_int[ ITYPCP ] == 2) lmin = 4;
     } //* fin zone parallele  *//



     //
     //
     //
     //* Declare memoire pour Newton ssiter_local
     //
     //
     //

    //taille du domaine a decouper (permet de prendre en compte la // de cond. limite)
    E_Int iv_lu =   ipt_param_int[ IJKV ] / size_ssdom[0]; E_Int jv_lu = ipt_param_int[ IJKV+1 ]/ size_ssdom[1];  E_Int kv_lu = ipt_param_int[ IJKV+2 ]  / size_ssdom[2];
    if(iv_lu== 0) iv_lu= 1;
    if(jv_lu== 0) jv_lu= 1;
    if(kv_lu== 0) kv_lu= 1;

    E_Int nisdom_lu  = iv_lu*jv_lu*kv_lu;
    E_Int ndimdx_dm  = 6*nisdom_lu*nssiter  // ind_dm
                      +            nssiter  // nisdom_residu
                      +            nssiter  // nidom_loc
                      +1                      ; // it_bloc

    ipt_param_int[ NDIMD_SDM ] = ndimdx_dm;
    ipt_param_int[ MXSSDOM_LU] = nisdom_lu;

    E_Int neq = 5;
    if (ipt_param_int[ IFLOW ] > 2) neq=6;

    //printf("iflow= %d %d %d %d \n", ipt_param_int[21], size_ssdom[0], size_ssdom[1], size_ssdom[2]);
    E_Int nitbuffer = 1;
    PyObject* rdm          = K_NUMPY::buildNumpyArray(nisdom_lu*neq*2*nssiter*nitbuffer, 1, 0, 1);
    PyObject* it_lu_ssdom  = K_NUMPY::buildNumpyArray(nisdom_lu*5            , 1, 1, 1);
    PyObject* ind_dm       = K_NUMPY::buildNumpyArray(ndimdx_dm              , 1, 1, 1);

    //E_Float* ipt_rdm          = K_NUMPY::getNumpyPtrF( rdm );
    E_Int*   ipt_ind_dm       = K_NUMPY::getNumpyPtrI( ind_dm );
    E_Int*   ipt_it_lu_ssdom  = K_NUMPY::getNumpyPtrI( it_lu_ssdom );

    
  //
  //Initialisation tableau  ind_dm(ndimdx_dm, nssiter)
  //
  for (E_Int nit = 0; nit < nssiter; nit++)  //boucle ss_iter rk3 ou Newton
   { 
    E_Int ndsdm = 0;
    for (E_Int k = 0; k < kv_lu; k++)   //!Boucle sous_domaine si Newton ssiter local
     { for (E_Int j = 0; j < jv_lu; j++)   
        { for (E_Int i = 0; i < iv_lu; i++)   
           {
             //E_Int ipos_ind_dm    = nit*nisdom_lu + ndsdm*6;
             E_Int ipos_ind_dm    = nit*nisdom_lu*6 + ndsdm*6;

             ipt_ind_dm[ ipos_ind_dm     ]= 1 + size_ssdom[0] *  i  ;
             ipt_ind_dm[ ipos_ind_dm + 1 ]=     size_ssdom[0] *(i+1);
             ipt_ind_dm[ ipos_ind_dm + 2 ]= 1 + size_ssdom[1] *  j  ;
             ipt_ind_dm[ ipos_ind_dm + 3 ]=     size_ssdom[1] *(j+1);
             ipt_ind_dm[ ipos_ind_dm + 4 ]= 1 + size_ssdom[2] *   k ;
             ipt_ind_dm[ ipos_ind_dm + 5 ]=     size_ssdom[2] *(k+1);
             if( i+1 == iv_lu ) ipt_ind_dm[ ipos_ind_dm + 1 ] = ipt_param_int[ IJKV  ];
             if( j+1 == jv_lu ) ipt_ind_dm[ ipos_ind_dm + 3 ] = ipt_param_int[ IJKV+1];
             if( k+1 == kv_lu ) ipt_ind_dm[ ipos_ind_dm + 5 ] = ipt_param_int[ IJKV+2];

             ndsdm++;
           }
        }
     }
   }
  //
  //Initialisation tableau  it_lu_ssdom( iv_lu*jv_lu*kv_lu=nisdom_lu  , nssiter)
  //
  for (E_Int ndsdm = 0; ndsdm < nisdom_lu; ndsdm++)   
  {
    ipt_it_lu_ssdom[    ndsdm                 ] = nssiter;  //it_lu_ssdom
    ipt_it_lu_ssdom[    ndsdm  + nisdom_lu    ] = nssiter;  //it_target_ssdom
    ipt_it_lu_ssdom[    ndsdm  + nisdom_lu*4  ] = ndsdm+1;  //no_lu
  }

  PyList_Append(metrics , ipti    );    Py_DECREF(ipti);
  PyList_Append(metrics , ipti0   );     
  PyList_Append(metrics , iptventi);     
  PyList_Append(metrics , ipti_df );     
  PyList_Append(metrics , rdm);         Py_DECREF(rdm);
  PyList_Append(metrics , ind_dm);      Py_DECREF(ind_dm );
  PyList_Append(metrics , it_lu_ssdom); Py_DECREF(it_lu_ssdom);
  PyList_Append(metrics , degener);     Py_DECREF(degener);

  if (lale>=1)   { Py_DECREF(ipti0); Py_DECREF(iptventi); }
  if (kfludom==3){ Py_DECREF(ipti_df); }

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
 
  return metrics;
}
