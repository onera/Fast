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
using namespace K_FLD;
using namespace std;

//=============================================================================
// Extrait les BC d'une zone d'un pyTree
// Il s'agit de donnees partagee avec python.
// Les donnees dont donc exactement les memes que celles du pyTree.
// IN: zone: la zone python
// OUT: nb_bc: nombre de CL sur la zone
// OUT: data_bc: data de la Cl
// OUT: ind_bc: fenetre de la CL
// OUT: type_bc: type de la CL (BCWall, BCFarfield,...)
// OUT: type_bcs: taille de la chaine type_bc
// Retourne 0: FAILED
// Sinon retourne le nombre de bcs de la zone.
//=============================================================================
E_Int K_FASTS::getbcfromzone(PyObject* zone ,
                             vector<E_Int>& data_pos ,
                             vector<E_Float*>& data_bc, 
                             vector<E_Int*>& ind_bc ,
                             vector<E_Int>& type_bc,
                             vector<PyArrayObject*>& hook) 
{
  PyObject* zonebc; E_Int nb_bc =0;

  zonebc = K_PYTREE::getNodeFromName1(zone, "ZoneBC");
  if (zonebc != NULL)
  {
    PyObject* list_bc = PyList_GetItem(zonebc, 2);
    nb_bc             = PyList_Size(list_bc);
    PyObject* bc; PyObject* node; char* str;
    
    //E_Int pos = 0;
    E_Int pos = data_bc.size();
    
    // boucle sur les noeuds contenus dans zone BC
    for (E_Int i = 0; i < nb_bc; i++)
    {
      bc   = PyList_GetItem(list_bc, i);
      node = PyList_GetItem(bc, 3);
      if (PyString_Check(node)) str  = PyString_AsString(node);  // type
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node)); 
#endif
      else str = NULL;
      if (K_STRING::cmp(str, "BC_t") == 0)
      {
        E_Int flag_data = 0;
        E_Int s, type=0;
        char* st = K_PYTREE::getValueS(bc, s, hook);
        // purely internal tag for BCs
        if      (K_STRING::cmp(st, s, "BCExtrapolate")           == 0) type = 0;
        else if (K_STRING::cmp(st, s, "BCFarfield")              == 0) type = 1;
        else if (K_STRING::cmp(st, s, "BCInflowSupersonic")      == 0) type = 2;
        else if (K_STRING::cmp(st, s, "BCWallViscous_transition")== 0) type =12;
        else if (K_STRING::cmp(st, s, "BCWallInviscid")          == 0) type = 3;
        else if (K_STRING::cmp(st, s, "BCSymmetryPlane")         == 0) type = 3;
        else if (K_STRING::cmp(st, s, "BCWall")                  == 0) type = 4;
        else if (K_STRING::cmp(st, s, "FamilySpecified")         == 0) type = 5;
        else if (K_STRING::cmp(st, s, "BCWallViscous")           == 0) type = 6;
        else if (K_STRING::cmp(st, s, "BCWallViscous_isot_fich") == 0) type = 7;
        else if (K_STRING::cmp(st, s, "BC_fich")                 == 0) type = 8;
        else if (K_STRING::cmp(st, s, "Nearmatch")               == 0) type = 9;
        else if (K_STRING::cmp(st, s, "BCOutflow")               == 0) type =10;
        else if (K_STRING::cmp(st, s, "BCInflow")                == 0) type =13;
        else if (K_STRING::cmp(st, s, "BCautoperiod")            == 0) type =11;
        else if (K_STRING::cmp(st, s, "BCDegenerateLine")        == 0) type = 0;
        else printf("Warning: getBCFromZone: unknown BC type %s.\n", st);

        type_bc.push_back(type);
        data_pos.push_back(pos);

       // printf("type_bc= %s %s %d \n", st, node[0] , type);
        
        PyObject* list_param = PyList_GetItem(bc, 2);
        E_Int nb_param       = PyList_Size(list_param);
        PyObject* param; 
        // boucle sur les noeuds parametre de la BC
        for (E_Int j = 0; j < nb_param; j++)
        {
          param = PyList_GetItem( list_param, j);
          node  = PyList_GetItem( param, 3);
          if (PyString_Check(node)) str = PyString_AsString(node);
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
#endif
          else str = NULL;
          if (K_STRING::cmp(str, "IndexRange_t") == 0)
          {
            E_Int s0, s1;
            E_Int* d = K_PYTREE::getValueAI(param, s0, s1, hook);
            ind_bc.push_back(d);
            if (s0*s1 != 6)
            {
              RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
              PyErr_SetString(PyExc_TypeError, "getbcfromzone: point range is invalid.");
              return 0; 
            }
          } 
          else if (K_STRING::cmp(str, "BCDataSet_t") == 0)
          {
            PyObject* list_dataset = PyList_GetItem(param, 2);
            E_Int nb_dataset       = PyList_Size(list_dataset);
            PyObject* dataset; 
            // boucle sur les noeuds data 
            for (E_Int k = 0; k < nb_dataset; k++)
            {
              dataset = PyList_GetItem( list_dataset, k);
              node    = PyList_GetItem( dataset, 3);
              if (PyString_Check(node)) str = PyString_AsString(node);
#if PY_VERSION_HEX >= 0x03000000
              else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
#endif
              else str = NULL;
              if (K_STRING::cmp(str, "BCData_t") == 0)
              {
                PyObject* list_datarray= PyList_GetItem(dataset, 2);
                E_Int nb_array         = PyList_Size(list_datarray);
                PyObject* datarray; 
                for (E_Int l = 0; l < nb_array; l++)
                {
                  datarray = PyList_GetItem(list_datarray, l);
                  node     = PyList_GetItem(datarray, 3);
                  if (PyString_Check(node)) str = PyString_AsString(node);
#if PY_VERSION_HEX >= 0x03000000
                  else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
#endif
                  else str = NULL;
                  if (K_STRING::cmp(str, "DataArray_t") == 0)
                  {
                    node = PyList_GetItem(datarray, 0); // var name
                    if (PyString_Check(node)) str = PyString_AsString(node);
#if PY_VERSION_HEX >= 0x03000000
                    else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
#endif
                    else str = NULL;
                    E_Float* f = K_PYTREE::getValueAF(datarray, hook);
                    data_bc.push_back(f);
                    
                    flag_data  = flag_data +1;
                  }
                  
                  
                }
                  
              }
            }
          } 
        } // boucle parametre
        
        if  (flag_data  == 0) { data_bc.push_back(NULL); pos  = pos + 1;}
        else                  {  pos  = pos + flag_data;}

      }//BC_t
    } //boucle noeud 
  }
  
  return nb_bc;
}
