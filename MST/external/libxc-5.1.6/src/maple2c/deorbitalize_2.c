v2rho2[0] = ked1_v2rho2[0]*mgga_vtau[0] + ked1_vrho[0]*ked1_vrho[0]*mgga_v2tau2[0] + 2*ked1_vrho[0]*mgga_v2rhotau[0] +
  mgga_v2rho2[0];
v2rhosigma[0] = ked1_v2rhosigma[0]*mgga_vtau[0] + ked1_vrho[0]*mgga_v2sigmatau[0] +
  ked1_vsigma[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) + mgga_v2rhosigma[0];
v2rholapl[0] = ked1_v2rholapl[0]*mgga_vtau[0] + ked1_vrho[0]*mgga_v2lapltau[0] +
  ked1_vlapl[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) + mgga_v2rholapl[0];
v2rhotau[0] = 0;
v2sigma2[0] = ked1_v2sigma2[0]*mgga_vtau[0] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v2tau2[0] +
  2*ked1_vsigma[0]*mgga_v2sigmatau[0] + mgga_v2sigma2[0];
v2sigmalapl[0] = ked1_v2sigmalapl[0]*mgga_vtau[0] + ked1_vsigma[0]*mgga_v2lapltau[0] +
  ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v2tau2[0] + mgga_v2sigmatau[0]) + mgga_v2sigmalapl[0];
v2sigmatau[0] = 0;
v2lapl2[0] = ked1_v2lapl2[0]*mgga_vtau[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v2tau2[0] +
  2*ked1_vlapl[0]*mgga_v2lapltau[0] + mgga_v2lapl2[0];
v2lapltau[0] = 0;
v2tau2[0] = 0;

if(func->nspin == XC_POLARIZED){
  v2rho2[1] = ked1_vrho[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) + ked2_vrho[0]*mgga_v2rhotau[1] +
    mgga_v2rho2[1];
  v2rho2[2] = ked2_v2rho2[0]*mgga_vtau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v2tau2[2] + 2*ked2_vrho[0]*mgga_v2rhotau[3]
    + mgga_v2rho2[2];
  v2rhosigma[1] = ked1_vrho[0]*mgga_v2sigmatau[2] + mgga_v2rhosigma[1];
  v2rhosigma[2] = ked1_vrho[0]*mgga_v2sigmatau[4] + ked2_vsigma[0]*(ked1_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[1]) +
    mgga_v2rhosigma[2];
  v2rhosigma[3] = ked2_vrho[0]*mgga_v2sigmatau[1] + ked1_vsigma[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    mgga_v2rhosigma[3];
  v2rhosigma[4] = ked2_vrho[0]*mgga_v2sigmatau[3] + mgga_v2rhosigma[4];
  v2rhosigma[5] = ked2_v2rhosigma[0]*mgga_vtau[1] + ked2_vrho[0]*mgga_v2sigmatau[5] +
    ked2_vsigma[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) + mgga_v2rhosigma[5];
  v2rholapl[1] = ked1_vrho[0]*mgga_v2lapltau[2] + ked2_vlapl[0]*(ked1_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[1]) +
    mgga_v2rholapl[1];
  v2rholapl[2] = ked2_vrho[0]*mgga_v2lapltau[1] + ked1_vlapl[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    mgga_v2rholapl[2];
  v2rholapl[3] = ked2_v2rholapl[0]*mgga_vtau[1] + ked2_vrho[0]*mgga_v2lapltau[3] +
    ked2_vlapl[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) + mgga_v2rholapl[3];
  v2rhotau[1] = 0;
  v2rhotau[2] = 0;
  v2rhotau[3] = 0;
  v2sigma2[1] = ked1_vsigma[0]*mgga_v2sigmatau[2] + mgga_v2sigma2[1];
  v2sigma2[2] = ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v2tau2[1] + mgga_v2sigmatau[4]) + ked2_vsigma[0]*mgga_v2sigmatau[1]
    + mgga_v2sigma2[2];
  v2sigma2[3] = mgga_v2sigma2[3];
  v2sigma2[4] = ked2_vsigma[0]*mgga_v2sigmatau[3] + mgga_v2sigma2[4];
  v2sigma2[5] = ked2_v2sigma2[0]*mgga_vtau[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v2tau2[2] +
    2*ked2_vsigma[0]*mgga_v2sigmatau[5] + mgga_v2sigma2[5];
  v2sigmalapl[1] = ked1_vsigma[0]*mgga_v2lapltau[2] + ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v2tau2[1] +
    mgga_v2sigmatau[1]) + mgga_v2sigmalapl[1];
  v2sigmalapl[2] = ked1_vlapl[0]*mgga_v2sigmatau[2] + mgga_v2sigmalapl[2];
  v2sigmalapl[3] = ked2_vlapl[0]*mgga_v2sigmatau[3] + mgga_v2sigmalapl[3];
  v2sigmalapl[4] = ked2_vsigma[0]*mgga_v2lapltau[1] + ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v2tau2[1] +
    mgga_v2sigmatau[4]) + mgga_v2sigmalapl[4];
  v2sigmalapl[5] = ked2_v2sigmalapl[0]*mgga_vtau[1] + ked2_vsigma[0]*mgga_v2lapltau[3] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v2tau2[2] + mgga_v2sigmatau[5]) + mgga_v2sigmalapl[5];
  v2sigmatau[1] = 0;
  v2sigmatau[2] = 0;
  v2sigmatau[3] = 0;
  v2sigmatau[4] = 0;
  v2sigmatau[5] = 0;
  v2lapl2[1] = ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v2tau2[1] + mgga_v2lapltau[2]) + ked2_vlapl[0]*mgga_v2lapltau[1] +
    mgga_v2lapl2[1];
  v2lapl2[2] = ked2_v2lapl2[0]*mgga_vtau[1] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v2tau2[2] +
    2*ked2_vlapl[0]*mgga_v2lapltau[3] + mgga_v2lapl2[2];
  v2lapltau[1] = 0;
  v2lapltau[2] = 0;
  v2lapltau[3] = 0;
  v2tau2[1] = 0;
  v2tau2[2] = 0;
}
