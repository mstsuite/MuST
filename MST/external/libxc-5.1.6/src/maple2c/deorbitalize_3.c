v3rho3[0] = ked1_v3rho3[0]*mgga_vtau[0] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[0] +
  3*ked1_v2rho2[0]*mgga_v2rhotau[0] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v3rhotau2[0] +
  3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v2tau2[0] + mgga_v3rho2tau[0]) + mgga_v3rho3[0];
v3rho2sigma[0] = ked1_v3rho2sigma[0]*mgga_vtau[0] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v2tau2[0] + mgga_v2sigmatau[0])
  + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3tau3[0] + mgga_v3sigmatau2[0]) +
  2*ked1_v2rhosigma[0]*mgga_v2rhotau[0] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[0] +
  ked1_vsigma[0]*mgga_v3rhotau2[0] + mgga_v3rhosigmatau[0]) + ked1_vsigma[0]*mgga_v3rho2tau[0] + mgga_v3rho2sigma[0];
v3rho2lapl[0] = ked1_v3rho2lapl[0]*mgga_vtau[0] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v2tau2[0] + mgga_v2lapltau[0]) +
  ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3tau3[0] + mgga_v3lapltau2[0]) + 2*ked1_v2rholapl[0]*mgga_v2rhotau[0]
  + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v2tau2[0] + ked1_vlapl[0]*mgga_v3rhotau2[0] + mgga_v3rholapltau[0]) +
  ked1_vlapl[0]*mgga_v3rho2tau[0] + mgga_v3rho2lapl[0];
v3rho2tau[0] = 0;
v3rhosigma2[0] = ked1_v3rhosigma2[0]*mgga_vtau[0] + 2*ked1_v2rhosigma[0]*mgga_v2sigmatau[0] +
  ked1_vrho[0]*mgga_v3sigma2tau[0] + ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) +
  ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vrho[0]*mgga_v3tau3[0] + mgga_v3rhotau2[0]) +
  2*ked1_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[0] + ked1_vrho[0]*mgga_v3sigmatau2[0] + mgga_v3rhosigmatau[0]) +
  mgga_v3rhosigma2[0];
v3rhosigmalapl[0] = ked1_v3rhosigmalapl[0]*mgga_vtau[0] + ked1_v2rhosigma[0]*(ked1_vlapl[0]*mgga_v2tau2[0] +
  mgga_v2lapltau[0]) + ked1_v2rholapl[0]*(ked1_vsigma[0]*mgga_v2tau2[0] + mgga_v2sigmatau[0]) +
  ked1_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[0] + ked1_vsigma[0]*mgga_v3lapltau2[0] +
  ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[0] + mgga_v3sigmatau2[0]) + mgga_v3sigmalapltau[0]) +
  ked1_v2sigmalapl[0]*mgga_v2rhotau[0] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3rhotau2[0] + mgga_v3rholapltau[0]) +
  ked1_vlapl[0]*mgga_v3rhosigmatau[0] + mgga_v3rhosigmalapl[0];
v3rhosigmatau[0] = 0;
v3rholapl2[0] = ked1_v3rholapl2[0]*mgga_vtau[0] + 2*ked1_v2rholapl[0]*mgga_v2lapltau[0] +
  ked1_vrho[0]*mgga_v3lapl2tau[0] + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) +
  ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vrho[0]*mgga_v3tau3[0] + mgga_v3rhotau2[0]) +
  2*ked1_vlapl[0]*(ked1_v2rholapl[0]*mgga_v2tau2[0] + ked1_vrho[0]*mgga_v3lapltau2[0] + mgga_v3rholapltau[0]) +
  mgga_v3rholapl2[0];
v3rholapltau[0] = 0;
v3rhotau2[0] = 0;
v3sigma3[0] = ked1_v3sigma3[0]*mgga_vtau[0] + ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[0] +
  3*ked1_v2sigma2[0]*mgga_v2sigmatau[0] + 3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3sigmatau2[0] +
  3*ked1_vsigma[0]*(ked1_v2sigma2[0]*mgga_v2tau2[0] + mgga_v3sigma2tau[0]) + mgga_v3sigma3[0];
v3sigma2lapl[0] = ked1_v3sigma2lapl[0]*mgga_vtau[0] + ked1_v2sigma2[0]*(ked1_vlapl[0]*mgga_v2tau2[0] +
  mgga_v2lapltau[0]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3tau3[0] + mgga_v3lapltau2[0]) +
  2*ked1_v2sigmalapl[0]*mgga_v2sigmatau[0] + 2*ked1_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[0] +
  ked1_vlapl[0]*mgga_v3sigmatau2[0] + mgga_v3sigmalapltau[0]) + ked1_vlapl[0]*mgga_v3sigma2tau[0] +
  mgga_v3sigma2lapl[0];
v3sigma2tau[0] = 0;
v3sigmalapl2[0] = ked1_v3sigmalapl2[0]*mgga_vtau[0] + 2*ked1_v2sigmalapl[0]*mgga_v2lapltau[0] +
  ked1_vsigma[0]*mgga_v3lapl2tau[0] + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v2tau2[0] + mgga_v2sigmatau[0]) +
  ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[0] + mgga_v3sigmatau2[0]) +
  2*ked1_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[0] + ked1_vsigma[0]*mgga_v3lapltau2[0] + mgga_v3sigmalapltau[0]) +
  mgga_v3sigmalapl2[0];
v3sigmalapltau[0] = 0;
v3sigmatau2[0] = 0;
v3lapl3[0] = ked1_v3lapl3[0]*mgga_vtau[0] + ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3tau3[0] +
  3*ked1_v2lapl2[0]*mgga_v2lapltau[0] + 3*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3lapltau2[0] +
  3*ked1_vlapl[0]*(ked1_v2lapl2[0]*mgga_v2tau2[0] + mgga_v3lapl2tau[0]) + mgga_v3lapl3[0];
v3lapl2tau[0] = 0;
v3lapltau2[0] = 0;
v3tau3[0] = 0;

if(func->nspin == XC_POLARIZED){
  v3rho3[1] = ked1_v2rho2[0]*mgga_v2rhotau[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v3rhotau2[3] +
    2*ked1_vrho[0]*mgga_v3rho2tau[2] + ked2_vrho[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) + mgga_v3rho3[1];
  v3rho3[2] = ked1_vrho[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] +
    2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) + ked2_v2rho2[0]*mgga_v2rhotau[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3rhotau2[2] + 2*ked2_vrho[0]*mgga_v3rho2tau[3] + mgga_v3rho3[2];
  v3rho3[3] = ked2_v3rho3[0]*mgga_vtau[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[3] +
    3*ked2_v2rho2[0]*mgga_v2rhotau[3] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v3rhotau2[5] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v2tau2[2] + mgga_v3rho2tau[5]) + mgga_v3rho3[3];
  v3rho2sigma[1] = ked1_v2rho2[0]*mgga_v2sigmatau[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v3sigmatau2[3] +
    2*ked1_vrho[0]*mgga_v3rhosigmatau[2] + mgga_v3rho2sigma[1];
  v3rho2sigma[2] = ked1_v2rho2[0]*mgga_v2sigmatau[4] + ked1_vrho[0]*ked1_vrho[0]*mgga_v3sigmatau2[6] +
    2*ked1_vrho[0]*mgga_v3rhosigmatau[4] + ked2_vsigma[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    mgga_v3rho2sigma[2];
  v3rho2sigma[3] = ked1_v2rhosigma[0]*mgga_v2rhotau[2] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3rhotau2[3] +
    mgga_v3rhosigmatau[6]) + ked2_vrho[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + ked1_vsigma[0]*mgga_v3rhotau2[1] +
    mgga_v3rhosigmatau[1]) + ked1_vsigma[0]*mgga_v3rho2tau[2] + mgga_v3rho2sigma[3];
  v3rho2sigma[4] = ked1_vrho[0]*(ked2_vrho[0]*mgga_v3sigmatau2[4] + mgga_v3rhosigmatau[8]) +
    ked2_vrho[0]*mgga_v3rhosigmatau[3] + mgga_v3rho2sigma[4];
  v3rho2sigma[5] = ked1_vrho[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[7]) + ked2_vsigma[0]*mgga_v3rhotau2[4] + mgga_v3rhosigmatau[10]) +
    ked2_v2rhosigma[0]*mgga_v2rhotau[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3rhotau2[2] + mgga_v3rhosigmatau[5]) +
    ked2_vsigma[0]*mgga_v3rho2tau[3] + mgga_v3rho2sigma[5];
  v3rho2sigma[6] = ked2_v2rho2[0]*mgga_v2sigmatau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3sigmatau2[2] +
    2*ked2_vrho[0]*mgga_v3rhosigmatau[7] + ked1_vsigma[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] + 2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) +
    mgga_v3rho2sigma[6];
  v3rho2sigma[7] = ked2_v2rho2[0]*mgga_v2sigmatau[3] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3sigmatau2[5] +
    2*ked2_vrho[0]*mgga_v3rhosigmatau[9] + mgga_v3rho2sigma[7];
  v3rho2sigma[8] = ked2_v3rho2sigma[0]*mgga_vtau[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v2tau2[2] +
    mgga_v2sigmatau[5]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3tau3[3] + mgga_v3sigmatau2[8]) +
    2*ked2_v2rhosigma[0]*mgga_v2rhotau[3] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[2] +
    ked2_vsigma[0]*mgga_v3rhotau2[5] + mgga_v3rhosigmatau[11]) + ked2_vsigma[0]*mgga_v3rho2tau[5] +
    mgga_v3rho2sigma[8];
  v3rho2lapl[1] = ked1_v2rho2[0]*mgga_v2lapltau[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v3lapltau2[3] +
    2*ked1_vrho[0]*mgga_v3rholapltau[2] + ked2_vlapl[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    mgga_v3rho2lapl[1];
  v3rho2lapl[2] = ked1_v2rholapl[0]*mgga_v2rhotau[2] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3rhotau2[3] +
    mgga_v3rholapltau[4]) + ked2_vrho[0]*(ked1_v2rholapl[0]*mgga_v2tau2[1] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3tau3[1]
    + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3rhotau2[1] + mgga_v3rholapltau[1]) + ked1_vlapl[0]*mgga_v3rho2tau[2] +
    mgga_v3rho2lapl[2];
  v3rho2lapl[3] = ked1_vrho[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3tau3[2] +
    mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3rhotau2[4] + mgga_v3rholapltau[6]) + ked2_v2rholapl[0]*mgga_v2rhotau[1]
    + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3rhotau2[2] + mgga_v3rholapltau[3]) + ked2_vlapl[0]*mgga_v3rho2tau[3] +
    mgga_v3rho2lapl[3];
  v3rho2lapl[4] = ked2_v2rho2[0]*mgga_v2lapltau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3lapltau2[2] +
    2*ked2_vrho[0]*mgga_v3rholapltau[5] + ked1_vlapl[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] + 2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) +
    mgga_v3rho2lapl[4];
  v3rho2lapl[5] = ked2_v3rho2lapl[0]*mgga_vtau[1] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v2tau2[2] + mgga_v2lapltau[3]) +
    ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3tau3[3] + mgga_v3lapltau2[5]) +
    2*ked2_v2rholapl[0]*mgga_v2rhotau[3] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*mgga_v3rhotau2[5] + mgga_v3rholapltau[7]) + ked2_vlapl[0]*mgga_v3rho2tau[5] + mgga_v3rho2lapl[5];
  v3rho2tau[1] = 0;
  v3rho2tau[2] = 0;
  v3rho2tau[3] = 0;
  v3rho2tau[4] = 0;
  v3rho2tau[5] = 0;
  v3rhosigma2[1] = ked1_v2rhosigma[0]*mgga_v2sigmatau[2] + ked1_vrho[0]*mgga_v3sigma2tau[2] +
    ked1_vsigma[0]*(ked1_vrho[0]*mgga_v3sigmatau2[3] + mgga_v3rhosigmatau[2]) + mgga_v3rhosigma2[1];
  v3rhosigma2[2] = ked1_v2rhosigma[0]*mgga_v2sigmatau[4] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[6] +
    mgga_v3sigma2tau[4]) + ked1_vsigma[0]*mgga_v3rhosigmatau[4] + ked2_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + ked1_vsigma[0]*mgga_v3rhotau2[1] +
    mgga_v3rhosigmatau[1]) + mgga_v3rhosigma2[2];
  v3rhosigma2[3] = ked1_vrho[0]*mgga_v3sigma2tau[6] + mgga_v3rhosigma2[3];
  v3rhosigma2[4] = ked1_vrho[0]*mgga_v3sigma2tau[8] + ked2_vsigma[0]*(ked1_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[3]) + mgga_v3rhosigma2[4];
  v3rhosigma2[5] = ked1_vrho[0]*mgga_v3sigma2tau[10] + ked2_v2sigma2[0]*(ked1_vrho[0]*mgga_v2tau2[1] +
    mgga_v2rhotau[1]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[2]) +
    2*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v3sigmatau2[7] + mgga_v3rhosigmatau[5]) + mgga_v3rhosigma2[5];
  v3rhosigma2[6] = ked2_vrho[0]*mgga_v3sigma2tau[1] + ked1_v2sigma2[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2])
    + ked1_vsigma[0]*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[3]) +
    2*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v3sigmatau2[1] + mgga_v3rhosigmatau[6]) + mgga_v3rhosigma2[6];
  v3rhosigma2[7] = ked2_vrho[0]*mgga_v3sigma2tau[3] + ked1_vsigma[0]*(ked2_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[8]) + mgga_v3rhosigma2[7];
  v3rhosigma2[8] = ked2_v2rhosigma[0]*mgga_v2sigmatau[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[2] +
    mgga_v3sigma2tau[5]) + ked1_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + ked2_vsigma[0]*mgga_v3rhotau2[4] +
    mgga_v3rhosigmatau[10]) + ked2_vsigma[0]*mgga_v3rhosigmatau[7] + mgga_v3rhosigma2[8];
  v3rhosigma2[9] = ked2_vrho[0]*mgga_v3sigma2tau[7] + mgga_v3rhosigma2[9];
  v3rhosigma2[10] = ked2_v2rhosigma[0]*mgga_v2sigmatau[3] + ked2_vrho[0]*mgga_v3sigma2tau[9] +
    ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3sigmatau2[5] + mgga_v3rhosigmatau[9]) + mgga_v3rhosigma2[10];
  v3rhosigma2[11] = ked2_v3rhosigma2[0]*mgga_vtau[1] + 2*ked2_v2rhosigma[0]*mgga_v2sigmatau[5] +
    ked2_vrho[0]*mgga_v3sigma2tau[11] + ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[3] + mgga_v3rhotau2[5]) +
    2*ked2_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[2] + ked2_vrho[0]*mgga_v3sigmatau2[8] + mgga_v3rhosigmatau[11]) +
    mgga_v3rhosigma2[11];
  v3rhosigmalapl[1] = ked1_v2rhosigma[0]*mgga_v2lapltau[2] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3lapltau2[3] +
    mgga_v3sigmalapltau[2]) + ked1_vsigma[0]*mgga_v3rholapltau[2] + ked2_vlapl[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + ked1_vsigma[0]*mgga_v3rhotau2[1] +
    mgga_v3rhosigmatau[1]) + mgga_v3rhosigmalapl[1];
  v3rhosigmalapl[2] = ked1_v2rholapl[0]*mgga_v2sigmatau[2] + ked1_vrho[0]*mgga_v3sigmalapltau[4] +
    ked1_vlapl[0]*(ked1_vrho[0]*mgga_v3sigmatau2[3] + mgga_v3rhosigmatau[2]) + mgga_v3rhosigmalapl[2];
  v3rhosigmalapl[3] = ked1_vrho[0]*mgga_v3sigmalapltau[6] + ked2_vlapl[0]*(ked1_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[3]) + mgga_v3rhosigmalapl[3];
  v3rhosigmalapl[4] = ked1_v2rholapl[0]*mgga_v2sigmatau[4] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[6] +
    mgga_v3sigmalapltau[8]) + ked2_vsigma[0]*(ked1_v2rholapl[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[1]) + ked1_vlapl[0]*mgga_v3rhosigmatau[4] + mgga_v3rhosigmalapl[4];
  v3rhosigmalapl[5] = ked1_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] + ked2_vsigma[0]*mgga_v3lapltau2[4] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + mgga_v3sigmalapltau[10]) +
    ked2_v2sigmalapl[0]*mgga_v2rhotau[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3rhotau2[2] + mgga_v3rholapltau[3]) +
    ked2_vlapl[0]*mgga_v3rhosigmatau[5] + mgga_v3rhosigmalapl[5];
  v3rhosigmalapl[6] = ked2_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[1] + ked1_vsigma[0]*mgga_v3lapltau2[1] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + mgga_v3sigmalapltau[1]) +
    ked1_v2sigmalapl[0]*mgga_v2rhotau[2] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3rhotau2[3] + mgga_v3rholapltau[4]) +
    ked1_vlapl[0]*mgga_v3rhosigmatau[6] + mgga_v3rhosigmalapl[6];
  v3rhosigmalapl[7] = ked2_v2rholapl[0]*mgga_v2sigmatau[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[2] +
    mgga_v3sigmalapltau[3]) + ked1_vsigma[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3rhotau2[4] +
    mgga_v3rholapltau[6]) + ked2_vlapl[0]*mgga_v3rhosigmatau[7] + mgga_v3rhosigmalapl[7];
  v3rhosigmalapl[8] = ked2_vrho[0]*mgga_v3sigmalapltau[5] + ked1_vlapl[0]*(ked2_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[8]) + mgga_v3rhosigmalapl[8];
  v3rhosigmalapl[9] = ked2_v2rholapl[0]*mgga_v2sigmatau[3] + ked2_vrho[0]*mgga_v3sigmalapltau[7] +
    ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3sigmatau2[5] + mgga_v3rhosigmatau[9]) + mgga_v3rhosigmalapl[9];
  v3rhosigmalapl[10] = ked2_v2rhosigma[0]*mgga_v2lapltau[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3lapltau2[2] +
    mgga_v3sigmalapltau[9]) + ked2_vsigma[0]*mgga_v3rholapltau[5] + ked1_vlapl[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + ked2_vsigma[0]*mgga_v3rhotau2[4] +
    mgga_v3rhosigmatau[10]) + mgga_v3rhosigmalapl[10];
  v3rhosigmalapl[11] = ked2_v3rhosigmalapl[0]*mgga_vtau[1] + ked2_v2rhosigma[0]*(ked2_vlapl[0]*mgga_v2tau2[2] +
    mgga_v2lapltau[3]) + ked2_v2rholapl[0]*(ked2_vsigma[0]*mgga_v2tau2[2] + mgga_v2sigmatau[5]) +
    ked2_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[2] + ked2_vsigma[0]*mgga_v3lapltau2[5] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[3] + mgga_v3sigmatau2[8]) + mgga_v3sigmalapltau[11]) +
    ked2_v2sigmalapl[0]*mgga_v2rhotau[3] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3rhotau2[5] + mgga_v3rholapltau[7]) +
    ked2_vlapl[0]*mgga_v3rhosigmatau[11] + mgga_v3rhosigmalapl[11];
  v3rhosigmatau[1] = 0;
  v3rhosigmatau[2] = 0;
  v3rhosigmatau[3] = 0;
  v3rhosigmatau[4] = 0;
  v3rhosigmatau[5] = 0;
  v3rhosigmatau[6] = 0;
  v3rhosigmatau[7] = 0;
  v3rhosigmatau[8] = 0;
  v3rhosigmatau[9] = 0;
  v3rhosigmatau[10] = 0;
  v3rhosigmatau[11] = 0;
  v3rholapl2[1] = ked1_v2rholapl[0]*mgga_v2lapltau[2] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3lapltau2[3] +
    mgga_v3lapl2tau[2]) + ked1_vlapl[0]*mgga_v3rholapltau[2] + ked2_vlapl[0]*(ked1_v2rholapl[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[1]) + mgga_v3rholapl2[1];
  v3rholapl2[2] = ked1_vrho[0]*mgga_v3lapl2tau[4] + ked2_v2lapl2[0]*(ked1_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[1]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[2]) +
    2*ked2_vlapl[0]*(ked1_vrho[0]*mgga_v3lapltau2[4] + mgga_v3rholapltau[3]) + mgga_v3rholapl2[2];
  v3rholapl2[3] = ked2_vrho[0]*mgga_v3lapl2tau[1] + ked1_v2lapl2[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[3]) +
    2*ked1_vlapl[0]*(ked2_vrho[0]*mgga_v3lapltau2[1] + mgga_v3rholapltau[4]) + mgga_v3rholapl2[3];
  v3rholapl2[4] = ked2_v2rholapl[0]*mgga_v2lapltau[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3lapltau2[2] +
    mgga_v3lapl2tau[3]) + ked1_vlapl[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3tau3[2]
    + mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3rhotau2[4] + mgga_v3rholapltau[6]) +
    ked2_vlapl[0]*mgga_v3rholapltau[5] + mgga_v3rholapl2[4];
  v3rholapl2[5] = ked2_v3rholapl2[0]*mgga_vtau[1] + 2*ked2_v2rholapl[0]*mgga_v2lapltau[3] +
    ked2_vrho[0]*mgga_v3lapl2tau[5] + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[3] + mgga_v3rhotau2[5]) +
    2*ked2_vlapl[0]*(ked2_v2rholapl[0]*mgga_v2tau2[2] + ked2_vrho[0]*mgga_v3lapltau2[5] + mgga_v3rholapltau[7]) +
    mgga_v3rholapl2[5];
  v3rholapltau[1] = 0;
  v3rholapltau[2] = 0;
  v3rholapltau[3] = 0;
  v3rholapltau[4] = 0;
  v3rholapltau[5] = 0;
  v3rholapltau[6] = 0;
  v3rholapltau[7] = 0;
  v3rhotau2[1] = 0;
  v3rhotau2[2] = 0;
  v3rhotau2[3] = 0;
  v3rhotau2[4] = 0;
  v3rhotau2[5] = 0;
  v3sigma3[1] = ked1_v2sigma2[0]*mgga_v2sigmatau[2] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3sigmatau2[3] +
    2*ked1_vsigma[0]*mgga_v3sigma2tau[2] + mgga_v3sigma3[1];
  v3sigma3[2] = ked1_v2sigma2[0]*mgga_v2sigmatau[4] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3sigmatau2[6] +
    2*ked1_vsigma[0]*mgga_v3sigma2tau[4] + ked2_vsigma[0]*(ked1_v2sigma2[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[1] + 2*ked1_vsigma[0]*mgga_v3sigmatau2[1] + mgga_v3sigma2tau[1]) +
    mgga_v3sigma3[2];
  v3sigma3[3] = ked1_vsigma[0]*mgga_v3sigma2tau[6] + mgga_v3sigma3[3];
  v3sigma3[4] = ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[4] + mgga_v3sigma2tau[8]) +
    ked2_vsigma[0]*mgga_v3sigma2tau[3] + mgga_v3sigma3[4];
  v3sigma3[5] = ked1_vsigma[0]*(ked2_v2sigma2[0]*mgga_v2tau2[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[2] +
    2*ked2_vsigma[0]*mgga_v3sigmatau2[7] + mgga_v3sigma2tau[10]) + ked2_v2sigma2[0]*mgga_v2sigmatau[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3sigmatau2[2] + 2*ked2_vsigma[0]*mgga_v3sigma2tau[5] + mgga_v3sigma3[5];
  v3sigma3[6] = mgga_v3sigma3[6];
  v3sigma3[7] = ked2_vsigma[0]*mgga_v3sigma2tau[7] + mgga_v3sigma3[7];
  v3sigma3[8] = ked2_v2sigma2[0]*mgga_v2sigmatau[3] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3sigmatau2[5] +
    2*ked2_vsigma[0]*mgga_v3sigma2tau[9] + mgga_v3sigma3[8];
  v3sigma3[9] = ked2_v3sigma3[0]*mgga_vtau[1] + ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[3] +
    3*ked2_v2sigma2[0]*mgga_v2sigmatau[5] + 3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3sigmatau2[8] +
    3*ked2_vsigma[0]*(ked2_v2sigma2[0]*mgga_v2tau2[2] + mgga_v3sigma2tau[11]) + mgga_v3sigma3[9];
  v3sigma2lapl[1] = ked1_v2sigma2[0]*mgga_v2lapltau[2] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3lapltau2[3] +
    2*ked1_vsigma[0]*mgga_v3sigmalapltau[2] + ked2_vlapl[0]*(ked1_v2sigma2[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[1] + 2*ked1_vsigma[0]*mgga_v3sigmatau2[1] + mgga_v3sigma2tau[1]) +
    mgga_v3sigma2lapl[1];
  v3sigma2lapl[2] = ked1_v2sigmalapl[0]*mgga_v2sigmatau[2] + ked1_vsigma[0]*mgga_v3sigmalapltau[4] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[3] + mgga_v3sigma2tau[2]) + mgga_v3sigma2lapl[2];
  v3sigma2lapl[3] = ked1_vsigma[0]*mgga_v3sigmalapltau[6] + ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[4] +
    mgga_v3sigma2tau[3]) + mgga_v3sigma2lapl[3];
  v3sigma2lapl[4] = ked1_v2sigmalapl[0]*mgga_v2sigmatau[4] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[6] +
    mgga_v3sigmalapltau[8]) + ked2_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3sigmatau2[1] +
    mgga_v3sigmalapltau[1]) + ked1_vlapl[0]*mgga_v3sigma2tau[4] + mgga_v3sigma2lapl[4];
  v3sigma2lapl[5] = ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3tau3[2] +
    mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3sigmatau2[7] + mgga_v3sigmalapltau[10]) +
    ked2_v2sigmalapl[0]*mgga_v2sigmatau[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[2] +
    mgga_v3sigmalapltau[3]) + ked2_vlapl[0]*mgga_v3sigma2tau[5] + mgga_v3sigma2lapl[5];
  v3sigma2lapl[6] = ked1_vlapl[0]*mgga_v3sigma2tau[6] + mgga_v3sigma2lapl[6];
  v3sigma2lapl[7] = ked2_vlapl[0]*mgga_v3sigma2tau[7] + mgga_v3sigma2lapl[7];
  v3sigma2lapl[8] = ked2_vsigma[0]*mgga_v3sigmalapltau[5] + ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[4] +
    mgga_v3sigma2tau[8]) + mgga_v3sigma2lapl[8];
  v3sigma2lapl[9] = ked2_v2sigmalapl[0]*mgga_v2sigmatau[3] + ked2_vsigma[0]*mgga_v3sigmalapltau[7] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[5] + mgga_v3sigma2tau[9]) + mgga_v3sigma2lapl[9];
  v3sigma2lapl[10] = ked2_v2sigma2[0]*mgga_v2lapltau[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3lapltau2[2] +
    2*ked2_vsigma[0]*mgga_v3sigmalapltau[9] + ked1_vlapl[0]*(ked2_v2sigma2[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[2] + 2*ked2_vsigma[0]*mgga_v3sigmatau2[7] + mgga_v3sigma2tau[10]) +
    mgga_v3sigma2lapl[10];
  v3sigma2lapl[11] = ked2_v3sigma2lapl[0]*mgga_vtau[1] + ked2_v2sigma2[0]*(ked2_vlapl[0]*mgga_v2tau2[2] +
    mgga_v2lapltau[3]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3tau3[3] + mgga_v3lapltau2[5]) +
    2*ked2_v2sigmalapl[0]*mgga_v2sigmatau[5] + 2*ked2_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*mgga_v3sigmatau2[8] + mgga_v3sigmalapltau[11]) + ked2_vlapl[0]*mgga_v3sigma2tau[11] +
    mgga_v3sigma2lapl[11];
  v3sigma2tau[1] = 0;
  v3sigma2tau[2] = 0;
  v3sigma2tau[3] = 0;
  v3sigma2tau[4] = 0;
  v3sigma2tau[5] = 0;
  v3sigma2tau[6] = 0;
  v3sigma2tau[7] = 0;
  v3sigma2tau[8] = 0;
  v3sigma2tau[9] = 0;
  v3sigma2tau[10] = 0;
  v3sigma2tau[11] = 0;
  v3sigmalapl2[1] = ked1_v2sigmalapl[0]*mgga_v2lapltau[2] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3lapltau2[3] +
    mgga_v3lapl2tau[2]) + ked1_vlapl[0]*mgga_v3sigmalapltau[2] + ked2_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3sigmatau2[1] +
    mgga_v3sigmalapltau[1]) + mgga_v3sigmalapl2[1];
  v3sigmalapl2[2] = ked1_vsigma[0]*mgga_v3lapl2tau[4] + ked2_v2lapl2[0]*(ked1_vsigma[0]*mgga_v2tau2[1] +
    mgga_v2sigmatau[1]) + ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[2]) +
    2*ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v3lapltau2[4] + mgga_v3sigmalapltau[3]) + mgga_v3sigmalapl2[2];
  v3sigmalapl2[3] = ked1_v2lapl2[0]*mgga_v2sigmatau[2] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3sigmatau2[3] +
    2*ked1_vlapl[0]*mgga_v3sigmalapltau[4] + mgga_v3sigmalapl2[3];
  v3sigmalapl2[4] = ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[4] + mgga_v3sigmalapltau[6]) +
    ked2_vlapl[0]*mgga_v3sigmalapltau[5] + mgga_v3sigmalapl2[4];
  v3sigmalapl2[5] = ked2_v2lapl2[0]*mgga_v2sigmatau[3] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3sigmatau2[5] +
    2*ked2_vlapl[0]*mgga_v3sigmalapltau[7] + mgga_v3sigmalapl2[5];
  v3sigmalapl2[6] = ked2_vsigma[0]*mgga_v3lapl2tau[1] + ked1_v2lapl2[0]*(ked2_vsigma[0]*mgga_v2tau2[1] +
    mgga_v2sigmatau[4]) + ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[6]) +
    2*ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v3lapltau2[1] + mgga_v3sigmalapltau[8]) + mgga_v3sigmalapl2[6];
  v3sigmalapl2[7] = ked2_v2sigmalapl[0]*mgga_v2lapltau[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3lapltau2[2] +
    mgga_v3lapl2tau[3]) + ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3sigmatau2[7] +
    mgga_v3sigmalapltau[10]) + ked2_vlapl[0]*mgga_v3sigmalapltau[9] + mgga_v3sigmalapl2[7];
  v3sigmalapl2[8] = ked2_v3sigmalapl2[0]*mgga_vtau[1] + 2*ked2_v2sigmalapl[0]*mgga_v2lapltau[3] +
    ked2_vsigma[0]*mgga_v3lapl2tau[5] + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v2tau2[2] + mgga_v2sigmatau[5]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[3] + mgga_v3sigmatau2[8]) +
    2*ked2_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[2] + ked2_vsigma[0]*mgga_v3lapltau2[5] + mgga_v3sigmalapltau[11])
    + mgga_v3sigmalapl2[8];
  v3sigmalapltau[1] = 0;
  v3sigmalapltau[2] = 0;
  v3sigmalapltau[3] = 0;
  v3sigmalapltau[4] = 0;
  v3sigmalapltau[5] = 0;
  v3sigmalapltau[6] = 0;
  v3sigmalapltau[7] = 0;
  v3sigmalapltau[8] = 0;
  v3sigmalapltau[9] = 0;
  v3sigmalapltau[10] = 0;
  v3sigmalapltau[11] = 0;
  v3sigmatau2[1] = 0;
  v3sigmatau2[2] = 0;
  v3sigmatau2[3] = 0;
  v3sigmatau2[4] = 0;
  v3sigmatau2[5] = 0;
  v3sigmatau2[6] = 0;
  v3sigmatau2[7] = 0;
  v3sigmatau2[8] = 0;
  v3lapl3[1] = ked1_v2lapl2[0]*mgga_v2lapltau[2] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3lapltau2[3] +
    2*ked1_vlapl[0]*mgga_v3lapl2tau[2] + ked2_vlapl[0]*(ked1_v2lapl2[0]*mgga_v2tau2[1] +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3tau3[1] + 2*ked1_vlapl[0]*mgga_v3lapltau2[1] + mgga_v3lapl2tau[1]) +
    mgga_v3lapl3[1];
  v3lapl3[2] = ked1_vlapl[0]*(ked2_v2lapl2[0]*mgga_v2tau2[1] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3tau3[2] +
    2*ked2_vlapl[0]*mgga_v3lapltau2[4] + mgga_v3lapl2tau[4]) + ked2_v2lapl2[0]*mgga_v2lapltau[1] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3lapltau2[2] + 2*ked2_vlapl[0]*mgga_v3lapl2tau[3] + mgga_v3lapl3[2];
  v3lapl3[3] = ked2_v3lapl3[0]*mgga_vtau[1] + ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3tau3[3] +
    3*ked2_v2lapl2[0]*mgga_v2lapltau[3] + 3*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3lapltau2[5] +
    3*ked2_vlapl[0]*(ked2_v2lapl2[0]*mgga_v2tau2[2] + mgga_v3lapl2tau[5]) + mgga_v3lapl3[3];
  v3lapl2tau[1] = 0;
  v3lapl2tau[2] = 0;
  v3lapl2tau[3] = 0;
  v3lapl2tau[4] = 0;
  v3lapl2tau[5] = 0;
  v3lapltau2[1] = 0;
  v3lapltau2[2] = 0;
  v3lapltau2[3] = 0;
  v3lapltau2[4] = 0;
  v3lapltau2[5] = 0;
  v3tau3[1] = 0;
  v3tau3[2] = 0;
  v3tau3[3] = 0;
}
