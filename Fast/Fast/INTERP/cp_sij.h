
E_Float dUdx, dUdy, dUdz, dVdx, dVdy, dVdz, dWdx, dWdy, dWdz, div_v;
E_Int   ip1,im1,jp1,jm1,kp1,km1, inci, incj, inck;

inci = 1;
incj = param_int[NoD][ NIJK];
inck = param_int[NoD][ NIJK]*param_int[NoD][ NIJK+1];

E_Float expo = 1./3.;
E_Float dx = pow(ipt_vol[NoD][0], expo);
//E_Float dx = ipt_vol[NoD][0]*3.;
E_Float dx_inv = 1./dx;

E_Float dz_inv = dx_inv;
if(bidim==1){dz_inv=0;}

if(param_int[NoD][ITYPZONE] != 2){printf("ERROR: calcul Sij pour couplage NS-LBM pas codé en non cartesien \n");}


for (E_Int noind = pt_deb; noind < pt_fin; noind++)
{

  indD0  = donorPts[noind];

  ip1 = indD0 + inci;
  im1 = indD0 - inci;
  jp1 = indD0 + incj;
  jm1 = indD0 - incj;
  kp1 = indD0 + inck;
  km1 = indD0 - inck;

  dUdx = 0.5*(vectOfDnrFields[1][ip1] - vectOfDnrFields[1][im1])*dx_inv;
  dUdy = 0.5*(vectOfDnrFields[1][jp1] - vectOfDnrFields[1][jm1])*dx_inv;
  dUdz = 0.5*(vectOfDnrFields[1][kp1] - vectOfDnrFields[1][km1])*dx_inv;

  dVdx = 0.5*(vectOfDnrFields[2][ip1] - vectOfDnrFields[2][im1])*dx_inv;
  dVdy = 0.5*(vectOfDnrFields[2][jp1] - vectOfDnrFields[2][jm1])*dx_inv;
  dVdz = 0.5*(vectOfDnrFields[2][kp1] - vectOfDnrFields[2][km1])*dx_inv;

  dWdx = 0.5*(vectOfDnrFields[3][ip1] - vectOfDnrFields[3][im1])*dz_inv;
  dWdy = 0.5*(vectOfDnrFields[3][jp1] - vectOfDnrFields[3][jm1])*dz_inv;
  dWdz = 0.5*(vectOfDnrFields[3][kp1] - vectOfDnrFields[3][km1])*dz_inv;

  div_v = dUdx + dVdy + dWdz;

  //printf("verif %.14f %d \n", dUdx - (2./3.)*div_v- vectOfDnrFields[ 5][indD0], indD0 );

  vectOfDnrFields[ 5][indD0] = dUdx - (2./3.)*div_v;
  vectOfDnrFields[ 8][indD0] = dVdy - (2./3.)*div_v;
  vectOfDnrFields[10][indD0] = dWdz - (2./3.)*div_v;
  vectOfDnrFields[ 6][indD0] = 0.5*(dUdy+dVdx);
  vectOfDnrFields[ 7][indD0] = 0.5*(dUdz+dWdx);
  vectOfDnrFields[ 9][indD0] = 0.5*(dVdz+dWdy);
}
