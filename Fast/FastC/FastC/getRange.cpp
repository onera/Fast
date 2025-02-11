/*   
    Copyright 2013-2025 Onera.

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
 
#include "FastC/fastc.h"
#include "FastC/param_solver.h"

using namespace std;
using namespace K_FLD;


//=============================================================================
PyObject* K_FASTC::PygetRange( PyObject* self, PyObject* args ) {
    PyObject* ind_bcArray;
    PyObject* param_intArray;
    E_Int     shift;
#if defined   E_DOUBLEINT
    if ( !PyArg_ParseTuple( args, "OOl", &ind_bcArray, &param_intArray, &shift ) ) return NULL;
#else
    if ( !PyArg_ParseTuple( args, "OOi", &ind_bcArray, &param_intArray, &shift ) ) return NULL;
#endif

    FldArrayI* ind_bc;
    FldArrayI* param_int;
    K_NUMPY::getFromNumpyArray( ind_bcArray, ind_bc, true );
    E_Int* ipt_ind_bc = ind_bc->begin( );
    K_NUMPY::getFromNumpyArray( param_intArray, param_int, true );
    E_Int* ipt_param_int = param_int->begin( );

    getRange( ipt_ind_bc, ipt_param_int + shift + 1, ipt_param_int );
    E_Int idir = getDir( ipt_param_int + IJKV, ipt_param_int + shift + 1 );

    ipt_param_int[shift] = idir;

    RELEASESHAREDN( ind_bcArray, ind_bc );
    RELEASESHAREDN( param_intArray, param_int );

    Py_INCREF( Py_None );
    return Py_None;
}
//=============================================================================
E_Int K_FASTC::getRange( E_Int* in_bc, E_Int* ind_fen, E_Int* param_int ) {
    // printf("adr  %d %d %d %d %d %d  \n", in_bc[0], in_bc[1],in_bc[2],in_bc[3],in_bc[4],in_bc[5] );
    if ( in_bc[0] == in_bc[1] ) {
        if ( in_bc[0] > 1 )
            ind_fen[0] = in_bc[0] - 1;
        else
            ind_fen[0] = in_bc[0];

        ind_fen[0] = std::max( E_Int(1), ind_fen[0] - 2 * param_int[NIJK + 3] );
        ind_fen[1] = ind_fen[0];

        ind_fen[2] = std::max( in_bc[2] - param_int[NIJK + 3], E_Int(1) );
        ind_fen[3] = std::min( in_bc[3] - param_int[NIJK + 3] - 1, param_int[IJKV + 1] );
        ind_fen[4] = std::max( in_bc[4] - param_int[NIJK + 4], E_Int(1) );
        ind_fen[5] = std::min( in_bc[5] - param_int[NIJK + 4] - 1, param_int[IJKV + 2] );
    } else if ( in_bc[2] == in_bc[3] ) {
        if ( in_bc[2] > 1 )
            ind_fen[2] = in_bc[2] - 1;
        else
            ind_fen[2] = in_bc[2];

        ind_fen[2] = std::max( E_Int(1), ind_fen[2] - 2 * param_int[NIJK + 3] );
        ind_fen[3] = ind_fen[2];

        ind_fen[0] = std::max( in_bc[0] - param_int[NIJK + 3], E_Int(1) );
        ind_fen[1] = std::min( in_bc[1] - param_int[NIJK + 3] - 1, param_int[IJKV] );
        ind_fen[4] = std::max( in_bc[4] - param_int[NIJK + 4], E_Int(1) );
        ind_fen[5] = std::min( in_bc[5] - param_int[NIJK + 4] - 1, param_int[IJKV + 2] );
    } else if ( in_bc[4] == in_bc[5] ) {
        if ( in_bc[4] > 1 )
            ind_fen[4] = in_bc[4] - 1;
        else
            ind_fen[4] = in_bc[4];

        ind_fen[4] = std::max( E_Int(1), ind_fen[4] - 2 * param_int[NIJK + 4] );
        ind_fen[5] = ind_fen[4];

        ind_fen[0] = std::max( in_bc[0] - param_int[NIJK + 3], E_Int(1) );
        ind_fen[1] = std::min( in_bc[1] - param_int[NIJK + 3] - 1, param_int[IJKV] );
        ind_fen[2] = std::max( in_bc[2] - param_int[NIJK + 3], E_Int(1) );
        ind_fen[3] = std::min( in_bc[3] - param_int[NIJK + 3] - 1, param_int[IJKV + 1] );
    }
     //printf("fen  %d %d %d %d %d %d  \n", ind_fen[0], ind_fen[1],ind_fen[2],ind_fen[3],ind_fen[4],ind_fen[5] );
     //printf("nijk  %d %d  \n", param_int[NIJK+3], param_int[NIJK+4] );
     //printf("ijkv  %d %d  %d\n", param_int[IJKV], param_int[IJKV+1], param_int[IJKV+2] );
    return 0;
}

//=============================================================================
E_Int K_FASTC::getDir( E_Int* ijkv, E_Int* ind_fen ) {
     //printf("%d %d %d %d %d %d %d %d %d \n", ind_fen[0],ind_fen[1],ind_fen[2],ind_fen[3],ind_fen[4],ind_fen[5],ijkv[0],ijkv[1],ijkv[2] );
    if ( ind_fen[0] == ind_fen[1] && ind_fen[0] == 1 ) return 1;
    if ( ind_fen[0] == ind_fen[1] && ind_fen[0] == ijkv[0] ) return 2;
    if ( ind_fen[2] == ind_fen[3] && ind_fen[2] == 1 ) return 3;
    if ( ind_fen[2] == ind_fen[3] && ind_fen[2] == ijkv[1] ) return 4;
    if ( ind_fen[4] == ind_fen[5] && ind_fen[4] == 1 ) return 5;
    if ( ind_fen[4] == ind_fen[5] && ind_fen[4] == ijkv[2] ) return 6;
    return 0;
}

