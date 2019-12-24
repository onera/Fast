/*    
    Copyright 2013-2020 Onera.

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
 
// Perfect gas laws

# include "fast.h"

//=============================================================================
// Constant gamma
//=============================================================================
E_Float K_FAST::gamma()
{
  return 1.4;
}
//=============================================================================
// Constant Prandtl 
//=============================================================================
E_Float K_FAST::prandtl()
{
  return 0.72;
}
//=============================================================================
// BetaS : constante de Sutherland = muS*(Ts+Cs)/(Ts*sqrt(Ts))
//=============================================================================
E_Float K_FAST::betaSuth(E_Float muSuth, E_Float CSuth, E_Float TSuth)
{
  return muSuth*(TSuth+CSuth)/(TSuth*sqrt(TSuth));
}
//=============================================================================
// Constant pressure specific heat 
//=============================================================================
E_Float K_FAST::Cp(E_Float Cv)
{
  return gamma()*Cv;
}
//=============================================================================
// Perfect Gas constant
//=============================================================================
E_Float K_FAST::Rgp(E_Float Cv)
{
  return (gamma()-1.)*Cv;
}
