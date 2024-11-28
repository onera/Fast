/*    
      Copyright 2013-2024 Onera.

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

# include "FastS/fastS.h"
# include "FastS/param_solver.h"
# include <string.h>
using namespace std;
using namespace K_FLD;

//=============================================================================
// Compute pour l'interface pyTree
//=============================================================================
PyObject* K_FASTS::display_ss_iteration(PyObject* self, PyObject* args)
{
  PyObject* zones; PyObject* metrics; PyObject* cvg_numpy;
  E_Int nitrun; E_Int lft; E_Int nssiter; E_Int iverb;
  
  if (!PYPARSETUPLE_(args, OO_ IIII_, &zones, &metrics, &nitrun, &nssiter, &lft, &iverb)) return NULL; 

  // tableau residu Newton de la zone
  E_Int neq_max = 7; E_Int niter_max = 100;
  FldArrayF rdm_glob(neq_max*niter_max*2); E_Float* iptrdm_glob = rdm_glob.begin(); 

  E_Int** ipt_ind_dm; E_Int** ipt_param_int; E_Float** iptrdm; E_Float** ipt_param_real;

  E_Int nidom = PyList_Size(zones);
  ipt_ind_dm      = new  E_Int*[2*nidom];
  ipt_param_int   = ipt_ind_dm  + nidom;

  iptrdm          = new  E_Float*[2*nidom];
  ipt_param_real  = iptrdm    + nidom;

  FldArrayF cvg(neq_max*4*nidom); E_Float* ipt_cvg = cvg.begin(); 
  vector<PyArrayObject*> hook;

  E_Float resLi[2]; resLi[0]=0.; resLi[1]=-1.e15; 
  for (E_Int nd = 0; nd < nidom; nd++)
  { 
    PyObject* zone = PyList_GetItem(zones, nd);
    char* name  = K_PYTREE::getNodeName(zone, hook); //nom de la zone
    E_Int size_name = strlen(name); // strlen retourne la longueur de la chaine sans le NULL final

    /* Get numerics from zone */
    PyObject*   numerics  = K_PYTREE::getNodeFromName1(zone    , ".Solver#ownData");
    PyObject*          o  = K_PYTREE::getNodeFromName1(numerics, "Parameter_int"); 
    ipt_param_int[nd]     = K_PYTREE::getValueAI(o, hook);
    o                     = K_PYTREE::getNodeFromName1(numerics, "Parameter_real"); 
    ipt_param_real[nd]    = K_PYTREE::getValueAF(o, hook);

    E_Int neq = ipt_param_int[nd][ NEQ ];
 
    E_Float* cvg_zone = ipt_cvg + 4*neq*nd;

    //E_Float* cvg_ptr; E_Int cvg_size; E_Int cvg_nfld;
    //K_NUMPY::getFromNumpyArray(cvg_numpy, cvg_ptr, cvg_size, cvg_nfld, true); 
 
    //printf(" cvg_size %d %d \n",cvg_size,cvg_nfld);

    if (lft >= 0)
    {
      PyObject* metric = PyList_GetItem(metrics, nd); // metric du domaine i

      // get metric
      iptrdm[ nd ]          =  K_NUMPY::getNumpyPtrF( PyList_GetItem(metric, METRIC_RDM ) );
      ipt_ind_dm[ nd ]      =  K_NUMPY::getNumpyPtrI( PyList_GetItem(metric, METRIC_INDM) );

      //
      //
      // on check si on a tout pour calculer sur la zone
      //
      //
      E_Int ijkv_lu[3];

      ijkv_lu[0] = K_FUNC::E_max( 1, ipt_param_int[nd][ IJKV    ]/ipt_param_int[nd][ SIZE_SSDOM   ]);
      ijkv_lu[1] = K_FUNC::E_max( 1, ipt_param_int[nd][ IJKV +1 ]/ipt_param_int[nd][ SIZE_SSDOM +1]);
      ijkv_lu[2] = K_FUNC::E_max( 1, ipt_param_int[nd][ IJKV +2 ]/ipt_param_int[nd][ SIZE_SSDOM +2]);

      //E_Int ndim_rdm= ijkv_lu[0]*ijkv_lu[1]*ijkv_lu[2];

      E_Int* ipt_it_bloc  = ipt_ind_dm[nd] + ipt_param_int[nd][ MXSSDOM_LU ]*6*nssiter + nssiter*2;    //it_bloc(nidom)
      E_Int  it_bloc_loc  = ipt_it_bloc[0];
      E_Float xpt_1       = 1./float( ipt_param_int[nd][ IJKV ]*ipt_param_int[nd][ IJKV +1 ]*ipt_param_int[nd][ IJKV +2 ]);

      PyObject* t2 = K_PYTREE::getNodeFromName1(numerics, "CFL_minmaxmoy");
      if (t2 != NULL) 
      { 
        E_Float* cfl = K_PYTREE::getValueAF(t2, hook);
        if (lft < 3) printf("CFL_MOY_MAX_MIN(Zone= %32.32s, nitrun=%d) %f    %f    %f \n", name, nitrun, cfl[2], cfl[0],cfl[1]); 
      }

      if (ipt_param_int[nd][ ITYPCP ] <= 1 || (ipt_param_int[nd][ ITYPCP ] == 2 && ipt_param_int[nd][ DTLOC ]==1))
      {
        if (it_bloc_loc*neq > niter_max*neq_max) printf("Display manometre: tableau trop petit");

        //assemblage norme Loo et L2
        for (E_Int ne= 0; ne < neq; ne++)
        { 
          for (E_Int i= 0; i < it_bloc_loc; i++)
          {
            E_Float rmax= 0.;
            E_Float rmoy= 0.;
            E_Int no_rdm;
            
            for (E_Int k_lu= 0; k_lu < ijkv_lu[2]; k_lu++) 
            {
              for (E_Int j_lu= 0; j_lu < ijkv_lu[1]; j_lu++) 
              {
                for (E_Int i_lu= 0; i_lu < ijkv_lu[0]; i_lu++) 
                {
                  no_rdm = i_lu + j_lu*ijkv_lu[0] + k_lu*ijkv_lu[0]*ijkv_lu[1];

                  E_Int ind1   = i + ne*nssiter +  no_rdm*neq*2*nssiter;
                  E_Int ind2   = i + (ne+neq)*nssiter +  no_rdm*neq*2*nssiter;

                  rmax   = max(rmax, iptrdm[nd][ind2] ); 
                  rmoy   = rmoy + iptrdm[nd][ind1];
                }
              }
            }
              
            E_Int ind1 =  it_bloc_loc*ne       + i;
            E_Int ind2 =  it_bloc_loc*(ne+neq) + i;

            rdm_glob[ind1] = sqrt(rmoy*xpt_1)/ipt_param_real[nd][ DTC ];
            rdm_glob[ind2] = rmax/ipt_param_real[nd][ DTC ];

            resLi[0] += rdm_glob[ind1];
            if (rdm_glob[ind2] > resLi[1]) resLi[1] =  rdm_glob[ind2];
          } //iteration
        } // neq

        dpssiter_(nitrun, neq,it_bloc_loc, nssiter, ipt_param_int[nd], lft, iverb, name, size_name, iptrdm_glob, cvg_zone );
      } // if zone implicit ou explicit+pdtloc
    } // if lft > 0

    /* Convergence History from zone */
    E_Float* ipt_Res_L2      ; E_Float* ipt_Res_oo      ; E_Int* ipt_Itnum ;
    E_Float* ipt_Res_L2_diff ; E_Float* ipt_Res_oo_diff ;
    E_Int LastRec;  E_Int nrec; E_Int nrec2;
    PyObject* zoneCH = K_PYTREE::getNodeFromName1(zone , "ZoneConvergenceHistory");
    if (zoneCH != NULL)
    {
      E_Int* ipt_LastRec = K_PYTREE::getValueAI(zoneCH, hook);
      LastRec = *ipt_LastRec;
      PyObject*  Itnum = K_PYTREE::getNodeFromName1(zoneCH, "IterationNumber");
      ipt_Itnum =  K_PYTREE::getValueAI(Itnum, nrec, nrec2, hook);
      PyObject*  Res_L2 = K_PYTREE::getNodeFromName1(zoneCH, "RSD_L2");
      ipt_Res_L2 =  K_PYTREE::getValueAF(Res_L2,hook) ;
      PyObject*  Res_oo = K_PYTREE::getNodeFromName1(zoneCH, "RSD_oo");
      ipt_Res_oo =  K_PYTREE::getValueAF(Res_oo, hook) ;

      PyObject*  Res_L2_diff = K_PYTREE::getNodeFromName1(zoneCH, "RSD_L2_diff");
      ipt_Res_L2_diff =  K_PYTREE::getValueAF(Res_L2_diff,hook) ;
      PyObject*  Res_oo_diff = K_PYTREE::getNodeFromName1(zoneCH, "RSD_oo_diff");
      ipt_Res_oo_diff =  K_PYTREE::getValueAF(Res_oo_diff, hook) ;

      //printf("display_ss_iteration : lft= %d nrec= %d\n",lft,nrec); fflush(stdout) ;
      conv2pytree_(cvg_zone, nitrun, neq, ipt_LastRec, name, size_name, lft, nrec, nd,ipt_Itnum, ipt_Res_L2, ipt_Res_oo,ipt_Res_L2_diff, ipt_Res_oo_diff);
    }  
  }// loop zone
  
  delete [] ipt_ind_dm;
  delete [] iptrdm;

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  PyObject* dico = PyDict_New();
  PyObject* tmp = Py_BuildValue("f", resLi[0]);
  PyDict_SetItemString(dico , "L2", tmp);
  tmp = Py_BuildValue("f", resLi[1]);
  PyDict_SetItemString(dico , "Loo", tmp);

  return dico; 
  
}
