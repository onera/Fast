// Les protos fortran
extern "C"
{  
  void cptimestepconv_(
    E_Int& cellNbTot, E_Float& cellDim,
    E_Float* ro, E_Float* rou, E_Float* rov, E_Float* row, E_Float* roE,
    E_Float* dtSteady,
    E_Float& gamma, E_Float& cfl);
  void cptimestepconvmotion_(
    E_Int& cellNbTot, 
    E_Float* sx, E_Float* sy, E_Float* sz,
    E_Float& cellDim,
    E_Float* ro, E_Float* rou, E_Float* rov, E_Float* row, E_Float* roE, 
    E_Float* dtSteady,
    E_Float& gamma, E_Float& cfl);
  void cptimesteptur_(
    E_Int &cellNbTot, E_Float* ro,
    E_Float& cellDim,
    E_Float* dtSteady, E_Float* timeStepConv,
    E_Float* mu, E_Float* muTurb,
    E_Float &gamma, E_Float& prandtl, E_Float& pranTurb, E_Float& cfl);
  void denconvec_(
    E_Int &cellt, E_Float* p,     
    E_Float* ro, E_Float* rou, E_Float* rov, E_Float* row, E_Float* roE,  
    E_Float* fcdx, E_Float* fcdy, E_Float* fcdz);
  void pressure_(
    E_Int& cellt,
    E_Float* ro, E_Float* rou, E_Float* rov, E_Float* row, E_Float* roE,
    E_Float& gam, E_Float* p);
  void cons2velocity_(
    E_Int& cellt, E_Float* ro, E_Float* rou, E_Float* rov, E_Float* row, 
    E_Float* u, E_Float* v, E_Float* w);
  void viscosity_(
    E_Int& cellt, E_Float& betas, E_Float& Cs, E_Float* temp, E_Float* mu);
  void heatcoef_(
    E_Int& cellt, E_Float& cp, E_Float& prandtl, E_Float* mu, E_Float* kappa);
  void temp_(
    E_Int& cellt, E_Float& cv, E_Float* ro, E_Float* u, E_Float* v, E_Float* w, E_Float* roE, E_Float* temp);
  void heatflux_(
    E_Int& cellt, E_Float* kappa, E_Float* gradxtemp,E_Float* gradytemp,E_Float* gradztemp,
    E_Float* qx, E_Float* qy, E_Float* qz);
  void visctensor_(
    E_Int& cellt, E_Float* mu,
    E_Float* gradxux, E_Float* gradxuy, E_Float* gradxuz,
    E_Float* gradyux, E_Float* gradyuy, E_Float* gradyuz,
    E_Float* gradzux, E_Float* gradzuy, E_Float* gradzuz,
    E_Float* tensxx, E_Float* tensxy, E_Float* tensxz,
    E_Float* tensyy, E_Float* tensyz, E_Float* tenszz);
  void densfluxdiff_(
    E_Int& cellt,
    E_Float* tgenxx, E_Float*  tgenxy, E_Float*  tgenxz,
    E_Float* tgenyy, E_Float*  tgenyz, E_Float*  tgenzz,
    E_Float* vax, E_Float*  vay, E_Float*  vaz, 
    E_Float* qgenx,  E_Float*  qgeny, E_Float*  qgenz,
    E_Float* fmomxdx, E_Float* fmomxdy, E_Float* fmomxdz,
    E_Float* fmomydx, E_Float* fmomydy, E_Float* fmomydz,
    E_Float* fmomzdx, E_Float* fmomzdy, E_Float* fmomzdz,
    E_Float* froedx,  E_Float* froedy,  E_Float* froedz);
  void setallbvaluesatf_(E_Int& length, E_Int& nfld, E_Float* x, 
                         E_Float& val, E_Int& ni, E_Int& nj, E_Int& nk, 
                         E_Int& border);
  void extrapallbvaluesf_(E_Int& length, E_Int& nfld, E_Float* x, 
                          E_Int& ni, E_Int& nj, E_Int& nk, 
                          E_Int& border);
}
