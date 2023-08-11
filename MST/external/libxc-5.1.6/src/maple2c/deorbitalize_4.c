v4rho4[0] = ked1_v4rho4[0]*mgga_vtau[0] + 3*ked1_v2rho2[0]*ked1_v2rho2[0]*mgga_v2tau2[0] +
  ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[0] + 4*ked1_v3rho3[0]*mgga_v2rhotau[0] +
  4*ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[0] +
  6*ked1_v2rho2[0]*(ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[0] + 2*ked1_vrho[0]*mgga_v3rhotau2[0] + mgga_v3rho2tau[0]) +
  6*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rho2tau2[0] + 4*ked1_vrho[0]*(ked1_v3rho3[0]*mgga_v2tau2[0] + mgga_v4rho3tau[0]) +
  mgga_v4rho4[0];
v4rho3sigma[0] = ked1_v4rho3sigma[0]*mgga_vtau[0] + ked1_vsigma[0]*ked1_v3rho3[0]*mgga_v2tau2[0] +
  ked1_v3rho3[0]*mgga_v2sigmatau[0] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4tau4[0] +
  mgga_v4sigmatau3[0]) + 3*ked1_v3rho2sigma[0]*mgga_v2rhotau[0] + 3*ked1_vsigma[0]*ked1_v2rho2[0]*mgga_v3rhotau2[0] +
  3*ked1_v2rho2[0]*mgga_v3rhosigmatau[0] + 3*ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4rhotau3[0] +
  mgga_v4rhosigmatau2[0]) + 3*ked1_v2rhosigma[0]*(ked1_v2rho2[0]*mgga_v2tau2[0] +
  ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[0] + 2*ked1_vrho[0]*mgga_v3rhotau2[0] + mgga_v3rho2tau[0]) +
  3*ked1_vrho[0]*(ked1_v3rho2sigma[0]*mgga_v2tau2[0] + ked1_v2rho2[0]*mgga_v3sigmatau2[0] +
  ked1_vsigma[0]*(ked1_v2rho2[0]*mgga_v3tau3[0] + mgga_v4rho2tau2[0]) + mgga_v4rho2sigmatau[0]) +
  ked1_vsigma[0]*mgga_v4rho3tau[0] + mgga_v4rho3sigma[0];
v4rho3lapl[0] = ked1_v4rho3lapl[0]*mgga_vtau[0] + ked1_vlapl[0]*ked1_v3rho3[0]*mgga_v2tau2[0] +
  ked1_v3rho3[0]*mgga_v2lapltau[0] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4tau4[0] +
  mgga_v4lapltau3[0]) + 3*ked1_v3rho2lapl[0]*mgga_v2rhotau[0] + 3*ked1_vlapl[0]*ked1_v2rho2[0]*mgga_v3rhotau2[0] +
  3*ked1_v2rho2[0]*mgga_v3rholapltau[0] + 3*ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4rhotau3[0] +
  mgga_v4rholapltau2[0]) + 3*ked1_v2rholapl[0]*(ked1_v2rho2[0]*mgga_v2tau2[0] +
  ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[0] + 2*ked1_vrho[0]*mgga_v3rhotau2[0] + mgga_v3rho2tau[0]) +
  3*ked1_vrho[0]*(ked1_v3rho2lapl[0]*mgga_v2tau2[0] + ked1_v2rho2[0]*mgga_v3lapltau2[0] +
  ked1_vlapl[0]*(ked1_v2rho2[0]*mgga_v3tau3[0] + mgga_v4rho2tau2[0]) + mgga_v4rho2lapltau[0]) +
  ked1_vlapl[0]*mgga_v4rho3tau[0] + mgga_v4rho3lapl[0];
v4rho3tau[0] = 0;
v4rho2sigma2[0] = ked1_v4rho2sigma2[0]*mgga_vtau[0] + 2*ked1_v2rhosigma[0]*ked1_v2rhosigma[0]*mgga_v2tau2[0] +
  ked1_v2sigma2[0]*ked1_v2rho2[0]*mgga_v2tau2[0] + 2*ked1_vsigma[0]*ked1_v3rho2sigma[0]*mgga_v2tau2[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*ked1_v2rho2[0]*mgga_v3tau3[0] + 2*ked1_v3rho2sigma[0]*mgga_v2sigmatau[0] +
  2*ked1_vsigma[0]*ked1_v2rho2[0]*mgga_v3sigmatau2[0] + ked1_v2rho2[0]*mgga_v3sigma2tau[0] +
  ked1_vrho[0]*ked1_vrho[0]*(ked1_v2sigma2[0]*mgga_v3tau3[0] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4tau4[0] +
  2*ked1_vsigma[0]*mgga_v4sigmatau3[0] + mgga_v4sigma2tau2[0]) + 2*ked1_v3rhosigma2[0]*mgga_v2rhotau[0] +
  4*ked1_v2rhosigma[0]*(ked1_vrho[0]*mgga_v3sigmatau2[0] + ked1_vsigma[0]*(ked1_vrho[0]*mgga_v3tau3[0] +
  mgga_v3rhotau2[0]) + mgga_v3rhosigmatau[0]) + 2*ked1_vrho[0]*(ked1_v3rhosigma2[0]*mgga_v2tau2[0] +
  ked1_v2sigma2[0]*mgga_v3rhotau2[0] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rhotau3[0] +
  2*ked1_vsigma[0]*mgga_v4rhosigmatau2[0] + mgga_v4rhosigma2tau[0]) + ked1_v2sigma2[0]*mgga_v3rho2tau[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rho2tau2[0] + 2*ked1_vsigma[0]*mgga_v4rho2sigmatau[0] + mgga_v4rho2sigma2[0];
v4rho2sigmalapl[0] = ked1_v4rho2sigmalapl[0]*mgga_vtau[0] + ked1_v2sigmalapl[0]*ked1_v2rho2[0]*mgga_v2tau2[0] +
  ked1_vsigma[0]*ked1_v3rho2lapl[0]*mgga_v2tau2[0] + ked1_vlapl[0]*ked1_v3rho2sigma[0]*mgga_v2tau2[0] +
  ked1_vlapl[0]*ked1_vsigma[0]*ked1_v2rho2[0]*mgga_v3tau3[0] + ked1_v3rho2sigma[0]*mgga_v2lapltau[0] +
  ked1_vsigma[0]*ked1_v2rho2[0]*mgga_v3lapltau2[0] + ked1_v3rho2lapl[0]*mgga_v2sigmatau[0] +
  ked1_vlapl[0]*ked1_v2rho2[0]*mgga_v3sigmatau2[0] + ked1_v2rho2[0]*mgga_v3sigmalapltau[0] +
  ked1_vrho[0]*ked1_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[0] + ked1_vsigma[0]*mgga_v4lapltau3[0] +
  ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[0] + mgga_v4sigmatau3[0]) + mgga_v4sigmalapltau2[0]) +
  2*ked1_v3rhosigmalapl[0]*mgga_v2rhotau[0] + 2*ked1_vlapl[0]*ked1_v2rhosigma[0]*mgga_v3rhotau2[0] +
  2*ked1_v2rhosigma[0]*mgga_v3rholapltau[0] + 2*ked1_v2rholapl[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[0] +
  ked1_vrho[0]*mgga_v3sigmatau2[0] + ked1_vsigma[0]*(ked1_vrho[0]*mgga_v3tau3[0] + mgga_v3rhotau2[0]) +
  mgga_v3rhosigmatau[0]) + 2*ked1_vrho[0]*(ked1_v3rhosigmalapl[0]*mgga_v2tau2[0] +
  ked1_v2rhosigma[0]*mgga_v3lapltau2[0] + ked1_v2sigmalapl[0]*mgga_v3rhotau2[0] + ked1_vsigma[0]*mgga_v4rholapltau2[0]
  + ked1_vlapl[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[0] + ked1_vsigma[0]*mgga_v4rhotau3[0] + mgga_v4rhosigmatau2[0]) +
  mgga_v4rhosigmalapltau[0]) + ked1_v2sigmalapl[0]*mgga_v3rho2tau[0] + ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4rho2tau2[0]
  + ked1_vsigma[0]*mgga_v4rho2lapltau[0] + ked1_vlapl[0]*mgga_v4rho2sigmatau[0] + mgga_v4rho2sigmalapl[0];
v4rho2sigmatau[0] = 0;
v4rho2lapl2[0] = ked1_v4rho2lapl2[0]*mgga_vtau[0] + 2*ked1_v2rholapl[0]*ked1_v2rholapl[0]*mgga_v2tau2[0] +
  ked1_v2lapl2[0]*ked1_v2rho2[0]*mgga_v2tau2[0] + 2*ked1_vlapl[0]*ked1_v3rho2lapl[0]*mgga_v2tau2[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_v2rho2[0]*mgga_v3tau3[0] + 2*ked1_v3rho2lapl[0]*mgga_v2lapltau[0] +
  2*ked1_vlapl[0]*ked1_v2rho2[0]*mgga_v3lapltau2[0] + ked1_v2rho2[0]*mgga_v3lapl2tau[0] +
  ked1_vrho[0]*ked1_vrho[0]*(ked1_v2lapl2[0]*mgga_v3tau3[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4tau4[0] +
  2*ked1_vlapl[0]*mgga_v4lapltau3[0] + mgga_v4lapl2tau2[0]) + 2*ked1_v3rholapl2[0]*mgga_v2rhotau[0] +
  4*ked1_v2rholapl[0]*(ked1_vrho[0]*mgga_v3lapltau2[0] + ked1_vlapl[0]*(ked1_vrho[0]*mgga_v3tau3[0] +
  mgga_v3rhotau2[0]) + mgga_v3rholapltau[0]) + 2*ked1_vrho[0]*(ked1_v3rholapl2[0]*mgga_v2tau2[0] +
  ked1_v2lapl2[0]*mgga_v3rhotau2[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rhotau3[0] +
  2*ked1_vlapl[0]*mgga_v4rholapltau2[0] + mgga_v4rholapl2tau[0]) + ked1_v2lapl2[0]*mgga_v3rho2tau[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rho2tau2[0] + 2*ked1_vlapl[0]*mgga_v4rho2lapltau[0] + mgga_v4rho2lapl2[0];
v4rho2lapltau[0] = 0;
v4rho2tau2[0] = 0;
v4rhosigma3[0] = ked1_v4rhosigma3[0]*mgga_vtau[0] + 3*ked1_vsigma[0]*ked1_v3rhosigma2[0]*mgga_v2tau2[0] +
  3*ked1_vsigma[0]*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3tau3[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4tau4[0] + 3*ked1_v3rhosigma2[0]*mgga_v2sigmatau[0] +
  6*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3sigmatau2[0] +
  3*ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigmatau3[0] + 3*ked1_v2rhosigma[0]*mgga_v3sigma2tau[0] +
  3*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigma2tau2[0] + ked1_vrho[0]*mgga_v4sigma3tau[0] +
  ked1_v3sigma3[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) +
  ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rhotau3[0] +
  3*ked1_v2sigma2[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[0] + ked1_vrho[0]*mgga_v3sigmatau2[0] +
  ked1_vsigma[0]*(ked1_vrho[0]*mgga_v3tau3[0] + mgga_v3rhotau2[0]) + mgga_v3rhosigmatau[0]) +
  3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rhosigmatau2[0] + 3*ked1_vsigma[0]*mgga_v4rhosigma2tau[0] +
  mgga_v4rhosigma3[0];
v4rhosigma2lapl[0] = ked1_v4rhosigma2lapl[0]*mgga_vtau[0] + 2*ked1_v2sigmalapl[0]*ked1_v2rhosigma[0]*mgga_v2tau2[0] +
  2*ked1_vsigma[0]*ked1_v3rhosigmalapl[0]*mgga_v2tau2[0] + ked1_vlapl[0]*ked1_v3rhosigma2[0]*mgga_v2tau2[0] +
  2*ked1_vsigma[0]*ked1_v2sigmalapl[0]*ked1_vrho[0]*mgga_v3tau3[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*ked1_v2rholapl[0]*mgga_v3tau3[0] +
  2*ked1_vlapl[0]*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3tau3[0] +
  ked1_vlapl[0]*ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4tau4[0] + ked1_v3rhosigma2[0]*mgga_v2lapltau[0] +
  2*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3lapltau2[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4lapltau3[0] + 2*ked1_v3rhosigmalapl[0]*mgga_v2sigmatau[0] +
  2*ked1_v2sigmalapl[0]*ked1_vrho[0]*mgga_v3sigmatau2[0] + 2*ked1_vsigma[0]*ked1_v2rholapl[0]*mgga_v3sigmatau2[0] +
  2*ked1_vlapl[0]*ked1_v2rhosigma[0]*mgga_v3sigmatau2[0] +
  2*ked1_vlapl[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigmatau3[0] + 2*ked1_v2rhosigma[0]*mgga_v3sigmalapltau[0] +
  2*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigmalapltau2[0] + ked1_v2rholapl[0]*mgga_v3sigma2tau[0] +
  ked1_vlapl[0]*ked1_vrho[0]*mgga_v4sigma2tau2[0] + ked1_vrho[0]*mgga_v4sigma2lapltau[0] +
  ked1_v3sigma2lapl[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) +
  2*ked1_vsigma[0]*ked1_v2sigmalapl[0]*mgga_v3rhotau2[0] +
  ked1_vlapl[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rhotau3[0] + ked1_v2sigma2[0]*(ked1_v2rholapl[0]*mgga_v2tau2[0] +
  ked1_vrho[0]*mgga_v3lapltau2[0] + ked1_vlapl[0]*(ked1_vrho[0]*mgga_v3tau3[0] + mgga_v3rhotau2[0]) +
  mgga_v3rholapltau[0]) + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rholapltau2[0] +
  2*ked1_v2sigmalapl[0]*mgga_v3rhosigmatau[0] + 2*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4rhosigmatau2[0] +
  2*ked1_vsigma[0]*mgga_v4rhosigmalapltau[0] + ked1_vlapl[0]*mgga_v4rhosigma2tau[0] + mgga_v4rhosigma2lapl[0];
v4rhosigma2tau[0] = 0;
v4rhosigmalapl2[0] = ked1_v4rhosigmalapl2[0]*mgga_vtau[0] + 2*ked1_v3rhosigmalapl[0]*(ked1_vlapl[0]*mgga_v2tau2[0] +
  mgga_v2lapltau[0]) + ked1_v2rhosigma[0]*(ked1_v2lapl2[0]*mgga_v2tau2[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3tau3[0]
  + 2*ked1_vlapl[0]*mgga_v3lapltau2[0] + mgga_v3lapl2tau[0]) + ked1_v3rholapl2[0]*(ked1_vsigma[0]*mgga_v2tau2[0] +
  mgga_v2sigmatau[0]) + 2*ked1_v2rholapl[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[0] + ked1_vsigma[0]*mgga_v3lapltau2[0] +
  ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[0] + mgga_v3sigmatau2[0]) + mgga_v3sigmalapltau[0]) +
  ked1_vrho[0]*(ked1_v3sigmalapl2[0]*mgga_v2tau2[0] + 2*ked1_v2sigmalapl[0]*mgga_v3lapltau2[0] +
  ked1_vsigma[0]*mgga_v4lapl2tau2[0] + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3tau3[0] + mgga_v3sigmatau2[0]) +
  ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[0] + mgga_v4sigmatau3[0]) +
  2*ked1_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[0] + ked1_vsigma[0]*mgga_v4lapltau3[0] + mgga_v4sigmalapltau2[0]) +
  mgga_v4sigmalapl2tau[0]) + ked1_v3sigmalapl2[0]*mgga_v2rhotau[0] +
  2*ked1_v2sigmalapl[0]*(ked1_vlapl[0]*mgga_v3rhotau2[0] + mgga_v3rholapltau[0]) +
  ked1_vsigma[0]*(ked1_v2lapl2[0]*mgga_v3rhotau2[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rhotau3[0] +
  2*ked1_vlapl[0]*mgga_v4rholapltau2[0] + mgga_v4rholapl2tau[0]) + ked1_v2lapl2[0]*mgga_v3rhosigmatau[0] +
  ked1_vlapl[0]*mgga_v4rhosigmalapltau[0] + ked1_vlapl[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[0] +
  mgga_v4rhosigmalapltau[0]) + mgga_v4rhosigmalapl2[0];
v4rhosigmalapltau[0] = 0;
v4rhosigmatau2[0] = 0;
v4rholapl3[0] = ked1_v4rholapl3[0]*mgga_vtau[0] + 3*ked1_vlapl[0]*ked1_v3rholapl2[0]*mgga_v2tau2[0] +
  3*ked1_vlapl[0]*ked1_vlapl[0]*ked1_v2rholapl[0]*mgga_v3tau3[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4tau4[0] + 3*ked1_v3rholapl2[0]*mgga_v2lapltau[0] +
  6*ked1_vlapl[0]*ked1_v2rholapl[0]*mgga_v3lapltau2[0] + 3*ked1_vlapl[0]*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4lapltau3[0]
  + 3*ked1_v2rholapl[0]*mgga_v3lapl2tau[0] + 3*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4lapl2tau2[0] +
  ked1_vrho[0]*mgga_v4lapl3tau[0] + ked1_v3lapl3[0]*(ked1_vrho[0]*mgga_v2tau2[0] + mgga_v2rhotau[0]) +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rhotau3[0] + 3*ked1_v2lapl2[0]*(ked1_v2rholapl[0]*mgga_v2tau2[0] +
  ked1_vrho[0]*mgga_v3lapltau2[0] + ked1_vlapl[0]*(ked1_vrho[0]*mgga_v3tau3[0] + mgga_v3rhotau2[0]) +
  mgga_v3rholapltau[0]) + 3*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rholapltau2[0] + 3*ked1_vlapl[0]*mgga_v4rholapl2tau[0] +
  mgga_v4rholapl3[0];
v4rholapl2tau[0] = 0;
v4rholapltau2[0] = 0;
v4rhotau3[0] = 0;
v4sigma4[0] = ked1_v4sigma4[0]*mgga_vtau[0] + 3*ked1_v2sigma2[0]*ked1_v2sigma2[0]*mgga_v2tau2[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4tau4[0] + 4*ked1_v3sigma3[0]*mgga_v2sigmatau[0] +
  4*ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[0] +
  6*ked1_v2sigma2[0]*(ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[0] + 2*ked1_vsigma[0]*mgga_v3sigmatau2[0] +
  mgga_v3sigma2tau[0]) + 6*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigma2tau2[0] +
  4*ked1_vsigma[0]*(ked1_v3sigma3[0]*mgga_v2tau2[0] + mgga_v4sigma3tau[0]) + mgga_v4sigma4[0];
v4sigma3lapl[0] = ked1_v4sigma3lapl[0]*mgga_vtau[0] + ked1_vlapl[0]*ked1_v3sigma3[0]*mgga_v2tau2[0] +
  ked1_v3sigma3[0]*mgga_v2lapltau[0] + ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[0] +
  mgga_v4lapltau3[0]) + 3*ked1_v3sigma2lapl[0]*mgga_v2sigmatau[0] +
  3*ked1_vlapl[0]*ked1_v2sigma2[0]*mgga_v3sigmatau2[0] + 3*ked1_v2sigma2[0]*mgga_v3sigmalapltau[0] +
  3*ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[0] + mgga_v4sigmalapltau2[0]) +
  3*ked1_v2sigmalapl[0]*(ked1_v2sigma2[0]*mgga_v2tau2[0] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[0] +
  2*ked1_vsigma[0]*mgga_v3sigmatau2[0] + mgga_v3sigma2tau[0]) + 3*ked1_vsigma[0]*(ked1_v3sigma2lapl[0]*mgga_v2tau2[0] +
  ked1_v2sigma2[0]*mgga_v3lapltau2[0] + ked1_vlapl[0]*(ked1_v2sigma2[0]*mgga_v3tau3[0] + mgga_v4sigma2tau2[0]) +
  mgga_v4sigma2lapltau[0]) + ked1_vlapl[0]*mgga_v4sigma3tau[0] + mgga_v4sigma3lapl[0];
v4sigma3tau[0] = 0;
v4sigma2lapl2[0] = ked1_v4sigma2lapl2[0]*mgga_vtau[0] + 2*ked1_v2sigmalapl[0]*ked1_v2sigmalapl[0]*mgga_v2tau2[0] +
  ked1_v2lapl2[0]*ked1_v2sigma2[0]*mgga_v2tau2[0] + 2*ked1_vlapl[0]*ked1_v3sigma2lapl[0]*mgga_v2tau2[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_v2sigma2[0]*mgga_v3tau3[0] + 2*ked1_v3sigma2lapl[0]*mgga_v2lapltau[0] +
  2*ked1_vlapl[0]*ked1_v2sigma2[0]*mgga_v3lapltau2[0] + ked1_v2sigma2[0]*mgga_v3lapl2tau[0] +
  ked1_vsigma[0]*ked1_vsigma[0]*(ked1_v2lapl2[0]*mgga_v3tau3[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4tau4[0] +
  2*ked1_vlapl[0]*mgga_v4lapltau3[0] + mgga_v4lapl2tau2[0]) + 2*ked1_v3sigmalapl2[0]*mgga_v2sigmatau[0] +
  4*ked1_v2sigmalapl[0]*(ked1_vsigma[0]*mgga_v3lapltau2[0] + ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[0] +
  mgga_v3sigmatau2[0]) + mgga_v3sigmalapltau[0]) + 2*ked1_vsigma[0]*(ked1_v3sigmalapl2[0]*mgga_v2tau2[0] +
  ked1_v2lapl2[0]*mgga_v3sigmatau2[0] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmatau3[0] +
  2*ked1_vlapl[0]*mgga_v4sigmalapltau2[0] + mgga_v4sigmalapl2tau[0]) + ked1_v2lapl2[0]*mgga_v3sigma2tau[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigma2tau2[0] + 2*ked1_vlapl[0]*mgga_v4sigma2lapltau[0] + mgga_v4sigma2lapl2[0];
v4sigma2lapltau[0] = 0;
v4sigma2tau2[0] = 0;
v4sigmalapl3[0] = ked1_v4sigmalapl3[0]*mgga_vtau[0] + 3*ked1_vlapl[0]*ked1_v3sigmalapl2[0]*mgga_v2tau2[0] +
  3*ked1_vlapl[0]*ked1_vlapl[0]*ked1_v2sigmalapl[0]*mgga_v3tau3[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4tau4[0] + 3*ked1_v3sigmalapl2[0]*mgga_v2lapltau[0] +
  6*ked1_vlapl[0]*ked1_v2sigmalapl[0]*mgga_v3lapltau2[0] +
  3*ked1_vlapl[0]*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4lapltau3[0] + 3*ked1_v2sigmalapl[0]*mgga_v3lapl2tau[0] +
  3*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4lapl2tau2[0] + ked1_vsigma[0]*mgga_v4lapl3tau[0] +
  ked1_v3lapl3[0]*(ked1_vsigma[0]*mgga_v2tau2[0] + mgga_v2sigmatau[0]) +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmatau3[0] + 3*ked1_v2lapl2[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[0]
  + ked1_vsigma[0]*mgga_v3lapltau2[0] + ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[0] + mgga_v3sigmatau2[0]) +
  mgga_v3sigmalapltau[0]) + 3*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmalapltau2[0] +
  3*ked1_vlapl[0]*mgga_v4sigmalapl2tau[0] + mgga_v4sigmalapl3[0];
v4sigmalapl2tau[0] = 0;
v4sigmalapltau2[0] = 0;
v4sigmatau3[0] = 0;
v4lapl4[0] = ked1_v4lapl4[0]*mgga_vtau[0] + 3*ked1_v2lapl2[0]*ked1_v2lapl2[0]*mgga_v2tau2[0] +
  ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4tau4[0] + 4*ked1_v3lapl3[0]*mgga_v2lapltau[0] +
  4*ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4lapltau3[0] +
  6*ked1_v2lapl2[0]*(ked1_vlapl[0]*ked1_vlapl[0]*mgga_v3tau3[0] + 2*ked1_vlapl[0]*mgga_v3lapltau2[0] +
  mgga_v3lapl2tau[0]) + 6*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4lapl2tau2[0] +
  4*ked1_vlapl[0]*(ked1_v3lapl3[0]*mgga_v2tau2[0] + mgga_v4lapl3tau[0]) + mgga_v4lapl4[0];
v4lapl3tau[0] = 0;
v4lapl2tau2[0] = 0;
v4lapltau3[0] = 0;
v4tau4[0] = 0;

if(func->nspin == XC_POLARIZED){
  v4rho4[1] = ked1_v3rho3[0]*mgga_v2rhotau[2] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[4] +
    3*ked1_v2rho2[0]*mgga_v3rho2tau[2] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rho2tau2[3] +
    3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3rhotau2[3] + mgga_v4rho3tau[2]) + ked2_vrho[0]*(ked1_v3rho3[0]*mgga_v2tau2[1]
    + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[1] + 3*ked1_v2rho2[0]*mgga_v3rhotau2[1] +
    3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[1] + 3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3tau3[1] + mgga_v4rho2tau2[1])
    + mgga_v4rho3tau[1]) + mgga_v4rho4[1];
  v4rho4[2] = ked1_v2rho2[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] +
    2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) + ked1_vrho[0]*ked1_vrho[0]*(ked2_v2rho2[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[2] + 2*ked2_vrho[0]*mgga_v4rhotau3[5] + mgga_v4rho2tau2[6]) +
    2*ked1_vrho[0]*(ked2_v2rho2[0]*mgga_v3rhotau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[2] +
    2*ked2_vrho[0]*mgga_v4rho2tau2[4] + mgga_v4rho3tau[4]) + ked2_v2rho2[0]*mgga_v3rho2tau[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4rho2tau2[2] + 2*ked2_vrho[0]*mgga_v4rho3tau[3] + mgga_v4rho4[2];
  v4rho4[3] = ked1_vrho[0]*(ked2_v3rho3[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[3] +
    3*ked2_v2rho2[0]*mgga_v3rhotau2[4] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[6] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3tau3[2] + mgga_v4rho2tau2[7]) + mgga_v4rho3tau[6]) +
    ked2_v3rho3[0]*mgga_v2rhotau[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[3] +
    3*ked2_v2rho2[0]*mgga_v3rho2tau[3] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rho2tau2[5] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3rhotau2[2] + mgga_v4rho3tau[5]) + mgga_v4rho4[3];
  v4rho4[4] = ked2_v4rho4[0]*mgga_vtau[1] + 3*ked2_v2rho2[0]*ked2_v2rho2[0]*mgga_v2tau2[2] +
    ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[4] + 4*ked2_v3rho3[0]*mgga_v2rhotau[3] +
    4*ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[7] +
    6*ked2_v2rho2[0]*(ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[3] + 2*ked2_vrho[0]*mgga_v3rhotau2[5] + mgga_v3rho2tau[5])
    + 6*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rho2tau2[8] + 4*ked2_vrho[0]*(ked2_v3rho3[0]*mgga_v2tau2[2] +
    mgga_v4rho3tau[7]) + mgga_v4rho4[4];
  v4rho3sigma[1] = ked1_v3rho3[0]*mgga_v2sigmatau[2] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[4] +
    3*ked1_v2rho2[0]*mgga_v3rhosigmatau[2] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhosigmatau2[3] +
    3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[3] + mgga_v4rho2sigmatau[2]) + mgga_v4rho3sigma[1];
  v4rho3sigma[2] = ked1_v3rho3[0]*mgga_v2sigmatau[4] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[8] +
    3*ked1_v2rho2[0]*mgga_v3rhosigmatau[4] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhosigmatau2[6] +
    3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[6] + mgga_v4rho2sigmatau[4]) +
    ked2_vsigma[0]*(ked1_v3rho3[0]*mgga_v2tau2[1] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[1] +
    3*ked1_v2rho2[0]*mgga_v3rhotau2[1] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[1] +
    3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3tau3[1] + mgga_v4rho2tau2[1]) + mgga_v4rho3tau[1]) + mgga_v4rho3sigma[2];
  v4rho3sigma[3] = ked1_v3rho2sigma[0]*mgga_v2rhotau[2] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3rhotau2[3] +
    mgga_v3rhosigmatau[6]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4rhotau3[4] + mgga_v4rhosigmatau2[9]) +
    2*ked1_v2rhosigma[0]*mgga_v3rho2tau[2] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3rhotau2[3] +
    ked1_vsigma[0]*mgga_v4rho2tau2[3] + mgga_v4rho2sigmatau[6]) + ked2_vrho[0]*(ked1_v3rho2sigma[0]*mgga_v2tau2[1] +
    ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) +
    ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) +
    2*ked1_v2rhosigma[0]*mgga_v3rhotau2[1] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[1] +
    ked1_vsigma[0]*mgga_v4rhotau3[1] + mgga_v4rhosigmatau2[1]) + ked1_vsigma[0]*mgga_v4rho2tau2[1] +
    mgga_v4rho2sigmatau[1]) + ked1_vsigma[0]*mgga_v4rho3tau[2] + mgga_v4rho3sigma[3];
  v4rho3sigma[4] = ked1_v2rho2[0]*mgga_v3rhosigmatau[8] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhosigmatau2[12] +
    2*ked1_vrho[0]*mgga_v4rho2sigmatau[8] + ked2_vrho[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[5] + 2*ked1_vrho[0]*mgga_v4rhosigmatau2[4] + mgga_v4rho2sigmatau[3]) +
    mgga_v4rho3sigma[4];
  v4rho3sigma[5] = ked2_vrho[0]*ked1_v2rho2[0]*mgga_v3sigmatau2[7] +
    ked1_vrho[0]*ked1_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[9] + ked1_v2rho2[0]*mgga_v3rhosigmatau[10] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhosigmatau2[15] + 2*ked1_vrho[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[7] +
    2*ked1_vrho[0]*mgga_v4rho2sigmatau[10] + ked2_v2rhosigma[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    ked2_vrho[0]*mgga_v4rho2sigmatau[5] + ked2_vsigma[0]*(ked1_v2rho2[0]*mgga_v3rhotau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[5] + 2*ked1_vrho[0]*mgga_v4rho2tau2[4] +
    ked2_vrho[0]*(ked1_v2rho2[0]*mgga_v3tau3[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[2] +
    2*ked1_vrho[0]*mgga_v4rhotau3[2] + mgga_v4rho2tau2[2]) + mgga_v4rho3tau[3]) + mgga_v4rho3sigma[5];
  v4rho3sigma[6] = ked1_vrho[0]*ked2_v2rho2[0]*mgga_v3sigmatau2[1] +
    ked1_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[2] + 2*ked1_vrho[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[10] +
    ked1_v2rhosigma[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] +
    2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) + ked1_vrho[0]*mgga_v4rho2sigmatau[12] +
    ked2_v2rho2[0]*mgga_v3rhosigmatau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[2] +
    2*ked2_vrho[0]*mgga_v4rho2sigmatau[7] + ked1_vsigma[0]*(ked1_vrho[0]*(ked2_v2rho2[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[2] + 2*ked2_vrho[0]*mgga_v4rhotau3[5] + mgga_v4rho2tau2[6]) +
    ked2_v2rho2[0]*mgga_v3rhotau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[2] + 2*ked2_vrho[0]*mgga_v4rho2tau2[4]
    + mgga_v4rho3tau[4]) + mgga_v4rho3sigma[6];
  v4rho3sigma[7] = ked1_vrho[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[4] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[6] +
    2*ked2_vrho[0]*mgga_v4rhosigmatau2[13] + mgga_v4rho2sigmatau[14]) + ked2_v2rho2[0]*mgga_v3rhosigmatau[3] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[5] + 2*ked2_vrho[0]*mgga_v4rho2sigmatau[9] + mgga_v4rho3sigma[7];
  v4rho3sigma[8] = ked1_vrho[0]*(ked2_v3rho2sigma[0]*mgga_v2tau2[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[7]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) +
    2*ked2_v2rhosigma[0]*mgga_v3rhotau2[4] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[2] +
    ked2_vsigma[0]*mgga_v4rhotau3[6] + mgga_v4rhosigmatau2[16]) + ked2_vsigma[0]*mgga_v4rho2tau2[7] +
    mgga_v4rho2sigmatau[16]) + ked2_v3rho2sigma[0]*mgga_v2rhotau[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3rhotau2[2]
    + mgga_v3rhosigmatau[5]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4rhotau3[3] + mgga_v4rhosigmatau2[8]) +
    2*ked2_v2rhosigma[0]*mgga_v3rho2tau[3] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3rhotau2[2] +
    ked2_vsigma[0]*mgga_v4rho2tau2[5] + mgga_v4rho2sigmatau[11]) + ked2_vsigma[0]*mgga_v4rho3tau[5] +
    mgga_v4rho3sigma[8];
  v4rho3sigma[9] = ked2_v3rho3[0]*mgga_v2sigmatau[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[3] +
    3*ked2_v2rho2[0]*mgga_v3rhosigmatau[7] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[11] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[2] + mgga_v4rho2sigmatau[13]) +
    ked1_vsigma[0]*(ked2_v3rho3[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[3] +
    3*ked2_v2rho2[0]*mgga_v3rhotau2[4] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[6] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3tau3[2] + mgga_v4rho2tau2[7]) + mgga_v4rho3tau[6]) + mgga_v4rho3sigma[9];
  v4rho3sigma[10] = ked2_v3rho3[0]*mgga_v2sigmatau[3] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[7] +
    3*ked2_v2rho2[0]*mgga_v3rhosigmatau[9] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[14] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[5] + mgga_v4rho2sigmatau[15]) + mgga_v4rho3sigma[10];
  v4rho3sigma[11] = ked2_v4rho3sigma[0]*mgga_vtau[1] + ked2_vsigma[0]*ked2_v3rho3[0]*mgga_v2tau2[2] +
    ked2_v3rho3[0]*mgga_v2sigmatau[5] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[4] +
    mgga_v4sigmatau3[11]) + 3*ked2_v3rho2sigma[0]*mgga_v2rhotau[3] + 3*ked2_vsigma[0]*ked2_v2rho2[0]*mgga_v3rhotau2[5]
    + 3*ked2_v2rho2[0]*mgga_v3rhosigmatau[11] + 3*ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4rhotau3[7] +
    mgga_v4rhosigmatau2[17]) + 3*ked2_v2rhosigma[0]*(ked2_v2rho2[0]*mgga_v2tau2[2] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[3] + 2*ked2_vrho[0]*mgga_v3rhotau2[5] + mgga_v3rho2tau[5]) +
    3*ked2_vrho[0]*(ked2_v3rho2sigma[0]*mgga_v2tau2[2] + ked2_v2rho2[0]*mgga_v3sigmatau2[8] +
    ked2_vsigma[0]*(ked2_v2rho2[0]*mgga_v3tau3[3] + mgga_v4rho2tau2[8]) + mgga_v4rho2sigmatau[17]) +
    ked2_vsigma[0]*mgga_v4rho3tau[7] + mgga_v4rho3sigma[11];
  v4rho3lapl[1] = ked1_v3rho3[0]*mgga_v2lapltau[2] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4lapltau3[4] +
    3*ked1_v2rho2[0]*mgga_v3rholapltau[2] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rholapltau2[3] +
    3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3lapltau2[3] + mgga_v4rho2lapltau[2]) +
    ked2_vlapl[0]*(ked1_v3rho3[0]*mgga_v2tau2[1] + ked1_vrho[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[1] +
    3*ked1_v2rho2[0]*mgga_v3rhotau2[1] + 3*ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[1] +
    3*ked1_vrho[0]*(ked1_v2rho2[0]*mgga_v3tau3[1] + mgga_v4rho2tau2[1]) + mgga_v4rho3tau[1]) + mgga_v4rho3lapl[1];
  v4rho3lapl[2] = ked1_v3rho2lapl[0]*mgga_v2rhotau[2] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3rhotau2[3] +
    mgga_v3rholapltau[4]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4rhotau3[4] + mgga_v4rholapltau2[6]) +
    2*ked1_v2rholapl[0]*mgga_v3rho2tau[2] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3rhotau2[3] +
    ked1_vlapl[0]*mgga_v4rho2tau2[3] + mgga_v4rho2lapltau[4]) + ked2_vrho[0]*(ked1_v3rho2lapl[0]*mgga_v2tau2[1] +
    ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) +
    ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) +
    2*ked1_v2rholapl[0]*mgga_v3rhotau2[1] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3tau3[1] +
    ked1_vlapl[0]*mgga_v4rhotau3[1] + mgga_v4rholapltau2[1]) + ked1_vlapl[0]*mgga_v4rho2tau2[1] +
    mgga_v4rho2lapltau[1]) + ked1_vlapl[0]*mgga_v4rho3tau[2] + mgga_v4rho3lapl[2];
  v4rho3lapl[3] = ked2_vrho[0]*ked1_v2rho2[0]*mgga_v3lapltau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*ked2_vrho[0]*mgga_v4lapltau3[5] + ked1_v2rho2[0]*mgga_v3rholapltau[6] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4rholapltau2[9] + 2*ked1_vrho[0]*ked2_vrho[0]*mgga_v4rholapltau2[4] +
    2*ked1_vrho[0]*mgga_v4rho2lapltau[6] + ked2_v2rholapl[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    ked2_vrho[0]*mgga_v4rho2lapltau[3] + ked2_vlapl[0]*(ked1_v2rho2[0]*mgga_v3rhotau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4rhotau3[5] + 2*ked1_vrho[0]*mgga_v4rho2tau2[4] +
    ked2_vrho[0]*(ked1_v2rho2[0]*mgga_v3tau3[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[2] +
    2*ked1_vrho[0]*mgga_v4rhotau3[2] + mgga_v4rho2tau2[2]) + mgga_v4rho3tau[3]) + mgga_v4rho3lapl[3];
  v4rho3lapl[4] = ked1_vrho[0]*ked2_v2rho2[0]*mgga_v3lapltau2[1] +
    ked1_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4lapltau3[2] + 2*ked1_vrho[0]*ked2_vrho[0]*mgga_v4rholapltau2[7] +
    ked1_v2rholapl[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] +
    2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) + ked1_vrho[0]*mgga_v4rho2lapltau[8] +
    ked2_v2rho2[0]*mgga_v3rholapltau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4rholapltau2[2] +
    2*ked2_vrho[0]*mgga_v4rho2lapltau[5] + ked1_vlapl[0]*(ked1_vrho[0]*(ked2_v2rho2[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[2] + 2*ked2_vrho[0]*mgga_v4rhotau3[5] + mgga_v4rho2tau2[6]) +
    ked2_v2rho2[0]*mgga_v3rhotau2[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[2] + 2*ked2_vrho[0]*mgga_v4rho2tau2[4]
    + mgga_v4rho3tau[4]) + mgga_v4rho3lapl[4];
  v4rho3lapl[5] = ked1_vrho[0]*(ked2_v3rho2lapl[0]*mgga_v2tau2[1] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3tau3[2] +
    mgga_v3lapltau2[4]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[3] + mgga_v4lapltau3[6]) +
    2*ked2_v2rholapl[0]*mgga_v3rhotau2[4] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3tau3[2] +
    ked2_vlapl[0]*mgga_v4rhotau3[6] + mgga_v4rholapltau2[10]) + ked2_vlapl[0]*mgga_v4rho2tau2[7] +
    mgga_v4rho2lapltau[10]) + ked2_v3rho2lapl[0]*mgga_v2rhotau[1] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3rhotau2[2] +
    mgga_v3rholapltau[3]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4rhotau3[3] + mgga_v4rholapltau2[5]) +
    2*ked2_v2rholapl[0]*mgga_v3rho2tau[3] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3rhotau2[2] +
    ked2_vlapl[0]*mgga_v4rho2tau2[5] + mgga_v4rho2lapltau[7]) + ked2_vlapl[0]*mgga_v4rho3tau[5] + mgga_v4rho3lapl[5];
  v4rho3lapl[6] = ked2_v3rho3[0]*mgga_v2lapltau[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4lapltau3[3] +
    3*ked2_v2rho2[0]*mgga_v3rholapltau[5] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rholapltau2[8] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3lapltau2[2] + mgga_v4rho2lapltau[9]) +
    ked1_vlapl[0]*(ked2_v3rho3[0]*mgga_v2tau2[1] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[3] +
    3*ked2_v2rho2[0]*mgga_v3rhotau2[4] + 3*ked2_vrho[0]*ked2_vrho[0]*mgga_v4rhotau3[6] +
    3*ked2_vrho[0]*(ked2_v2rho2[0]*mgga_v3tau3[2] + mgga_v4rho2tau2[7]) + mgga_v4rho3tau[6]) + mgga_v4rho3lapl[6];
  v4rho3lapl[7] = ked2_v4rho3lapl[0]*mgga_vtau[1] + ked2_vlapl[0]*ked2_v3rho3[0]*mgga_v2tau2[2] +
    ked2_v3rho3[0]*mgga_v2lapltau[3] + ked2_vrho[0]*ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[4] +
    mgga_v4lapltau3[7]) + 3*ked2_v3rho2lapl[0]*mgga_v2rhotau[3] + 3*ked2_vlapl[0]*ked2_v2rho2[0]*mgga_v3rhotau2[5] +
    3*ked2_v2rho2[0]*mgga_v3rholapltau[7] + 3*ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4rhotau3[7] +
    mgga_v4rholapltau2[11]) + 3*ked2_v2rholapl[0]*(ked2_v2rho2[0]*mgga_v2tau2[2] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[3] + 2*ked2_vrho[0]*mgga_v3rhotau2[5] + mgga_v3rho2tau[5]) +
    3*ked2_vrho[0]*(ked2_v3rho2lapl[0]*mgga_v2tau2[2] + ked2_v2rho2[0]*mgga_v3lapltau2[5] +
    ked2_vlapl[0]*(ked2_v2rho2[0]*mgga_v3tau3[3] + mgga_v4rho2tau2[8]) + mgga_v4rho2lapltau[11]) +
    ked2_vlapl[0]*mgga_v4rho3tau[7] + mgga_v4rho3lapl[7];
  v4rho3tau[1] = 0;
  v4rho3tau[2] = 0;
  v4rho3tau[3] = 0;
  v4rho3tau[4] = 0;
  v4rho3tau[5] = 0;
  v4rho3tau[6] = 0;
  v4rho3tau[7] = 0;
  v4rho2sigma2[1] = ked1_v3rho2sigma[0]*mgga_v2sigmatau[2] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[3] +
    mgga_v3sigma2tau[2]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[4] + mgga_v4sigma2tau2[3]) +
    2*ked1_v2rhosigma[0]*mgga_v3rhosigmatau[2] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[3] +
    ked1_vsigma[0]*mgga_v4rhosigmatau2[3] + mgga_v4rhosigma2tau[2]) + ked1_vsigma[0]*mgga_v4rho2sigmatau[2] +
    mgga_v4rho2sigma2[1];
  v4rho2sigma2[2] = ked1_v3rho2sigma[0]*mgga_v2sigmatau[4] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[6] +
    mgga_v3sigma2tau[4]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[8] + mgga_v4sigma2tau2[6]) +
    2*ked1_v2rhosigma[0]*mgga_v3rhosigmatau[4] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[6] +
    ked1_vsigma[0]*mgga_v4rhosigmatau2[6] + mgga_v4rhosigma2tau[4]) + ked1_vsigma[0]*mgga_v4rho2sigmatau[4] +
    ked2_vsigma[0]*(ked1_v3rho2sigma[0]*mgga_v2tau2[1] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3tau3[1] +
    mgga_v3sigmatau2[1]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) +
    2*ked1_v2rhosigma[0]*mgga_v3rhotau2[1] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[1] +
    ked1_vsigma[0]*mgga_v4rhotau3[1] + mgga_v4rhosigmatau2[1]) + ked1_vsigma[0]*mgga_v4rho2tau2[1] +
    mgga_v4rho2sigmatau[1]) + mgga_v4rho2sigma2[2];
  v4rho2sigma2[3] = ked1_v2rho2[0]*mgga_v3sigma2tau[6] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigma2tau2[9] +
    2*ked1_vrho[0]*mgga_v4rhosigma2tau[6] + mgga_v4rho2sigma2[3];
  v4rho2sigma2[4] = ked1_v2rho2[0]*mgga_v3sigma2tau[8] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigma2tau2[12] +
    2*ked1_vrho[0]*mgga_v4rhosigma2tau[8] + ked2_vsigma[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[5] + 2*ked1_vrho[0]*mgga_v4rhosigmatau2[4] + mgga_v4rho2sigmatau[3]) +
    mgga_v4rho2sigma2[4];
  v4rho2sigma2[5] = ked1_v2rho2[0]*mgga_v3sigma2tau[10] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigma2tau2[15] +
    2*ked1_vrho[0]*mgga_v4rhosigma2tau[10] + ked2_v2sigma2[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked1_v2rho2[0]*mgga_v3tau3[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[2] +
    2*ked1_vrho[0]*mgga_v4rhotau3[2] + mgga_v4rho2tau2[2]) + 2*ked2_vsigma[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[7] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[9] + 2*ked1_vrho[0]*mgga_v4rhosigmatau2[7] + mgga_v4rho2sigmatau[5]) +
    mgga_v4rho2sigma2[5];
  v4rho2sigma2[6] = ked1_v3rhosigma2[0]*mgga_v2rhotau[2] + 2*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3rhotau2[3] +
    ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4rhotau3[4] + 2*ked1_v2rhosigma[0]*mgga_v3rhosigmatau[6] +
    2*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4rhosigmatau2[9] + ked1_vrho[0]*mgga_v4rhosigma2tau[12] +
    ked2_vrho[0]*(ked1_v3rhosigma2[0]*mgga_v2tau2[1] + 2*ked1_v2rhosigma[0]*mgga_v3sigmatau2[1] +
    ked1_vrho[0]*mgga_v4sigma2tau2[1] + ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[1]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[1]) +
    2*ked1_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[1] + ked1_vrho[0]*mgga_v4sigmatau3[1] + mgga_v4rhosigmatau2[1]) +
    mgga_v4rhosigma2tau[1]) + ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3rhotau2[3] + mgga_v3rho2tau[2]) +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rho2tau2[3] + 2*ked1_vsigma[0]*mgga_v4rho2sigmatau[6] + mgga_v4rho2sigma2[6];
  v4rho2sigma2[7] = ked1_v2rhosigma[0]*mgga_v3rhosigmatau[8] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4rhosigmatau2[12] +
    mgga_v4rhosigma2tau[14]) + ked2_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[4] +
    ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[5] + mgga_v4sigma2tau2[4]) + ked1_vsigma[0]*mgga_v4rhosigmatau2[4] +
    mgga_v4rhosigma2tau[3]) + ked1_vsigma[0]*mgga_v4rho2sigmatau[8] + mgga_v4rho2sigma2[7];
  v4rho2sigma2[8] = ked1_v2rhosigma[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] + ked2_vrho[0]*mgga_v3sigmatau2[7] +
    ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) + mgga_v3rhosigmatau[10]) +
    ked1_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[2] +
    mgga_v4sigma2tau2[7]) + ked1_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[2] + mgga_v4sigmatau3[9]) + ked2_vsigma[0]*mgga_v4rhotau3[5] +
    mgga_v4rhosigmatau2[15]) + ked2_vsigma[0]*mgga_v4rhosigmatau2[10] + mgga_v4rhosigma2tau[16]) +
    ked2_v2rhosigma[0]*(ked1_vsigma[0]*mgga_v3rhotau2[1] + mgga_v3rhosigmatau[1]) +
    ked2_vrho[0]*(ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v4rhotau3[2] + mgga_v4rhosigmatau2[7]) +
    ked2_vsigma[0]*mgga_v4rhosigmatau2[2] + mgga_v4rhosigma2tau[5]) + ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v4rho2tau2[4]
    + mgga_v4rho2sigmatau[10]) + ked2_vsigma[0]*mgga_v4rho2sigmatau[7] + mgga_v4rho2sigma2[8];
  v4rho2sigma2[9] = ked1_vrho[0]*(ked2_vrho[0]*mgga_v4sigma2tau2[10] + mgga_v4rhosigma2tau[18]) +
    ked2_vrho[0]*mgga_v4rhosigma2tau[7] + mgga_v4rho2sigma2[9];
  v4rho2sigma2[10] = ked1_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[6] + mgga_v4sigma2tau2[13]) + ked2_vsigma[0]*mgga_v4rhosigmatau2[13]
    + mgga_v4rhosigma2tau[20]) + ked2_v2rhosigma[0]*mgga_v3rhosigmatau[3] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4rhosigmatau2[5] + mgga_v4rhosigma2tau[9]) +
    ked2_vsigma[0]*mgga_v4rho2sigmatau[9] + mgga_v4rho2sigma2[10];
  v4rho2sigma2[11] = ked1_vrho[0]*(ked2_v3rhosigma2[0]*mgga_v2tau2[1] + 2*ked2_v2rhosigma[0]*mgga_v3sigmatau2[7] +
    ked2_vrho[0]*mgga_v4sigma2tau2[16] + ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[6]) +
    2*ked2_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[2] + ked2_vrho[0]*mgga_v4sigmatau3[10] + mgga_v4rhosigmatau2[16])
    + mgga_v4rhosigma2tau[22]) + ked2_v3rhosigma2[0]*mgga_v2rhotau[1] +
    2*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3rhotau2[2] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4rhotau3[3] + 2*ked2_v2rhosigma[0]*mgga_v3rhosigmatau[5] +
    2*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4rhosigmatau2[8] + ked2_vrho[0]*mgga_v4rhosigma2tau[11] +
    ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3rhotau2[2] + mgga_v3rho2tau[3]) +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rho2tau2[5] + 2*ked2_vsigma[0]*mgga_v4rho2sigmatau[11] +
    mgga_v4rho2sigma2[11];
  v4rho2sigma2[12] = ked2_v2rho2[0]*mgga_v3sigma2tau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigma2tau2[2] +
    2*ked2_vrho[0]*mgga_v4rhosigma2tau[13] + ked1_v2sigma2[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] + 2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked2_v2rho2[0]*mgga_v3tau3[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[2] +
    2*ked2_vrho[0]*mgga_v4rhotau3[5] + mgga_v4rho2tau2[6]) + 2*ked1_vsigma[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[2] + 2*ked2_vrho[0]*mgga_v4rhosigmatau2[10] + mgga_v4rho2sigmatau[12]) +
    mgga_v4rho2sigma2[12];
  v4rho2sigma2[13] = ked2_v2rho2[0]*mgga_v3sigma2tau[3] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigma2tau2[5] +
    2*ked2_vrho[0]*mgga_v4rhosigma2tau[15] + ked1_vsigma[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[6] + 2*ked2_vrho[0]*mgga_v4rhosigmatau2[13] + mgga_v4rho2sigmatau[14]) +
    mgga_v4rho2sigma2[13];
  v4rho2sigma2[14] = ked2_v3rho2sigma[0]*mgga_v2sigmatau[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[2] +
    mgga_v3sigma2tau[5]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[3] + mgga_v4sigma2tau2[8]) +
    2*ked2_v2rhosigma[0]*mgga_v3rhosigmatau[7] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[2] +
    ked2_vsigma[0]*mgga_v4rhosigmatau2[11] + mgga_v4rhosigma2tau[17]) +
    ked1_vsigma[0]*(ked2_v3rho2sigma[0]*mgga_v2tau2[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[7]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) +
    2*ked2_v2rhosigma[0]*mgga_v3rhotau2[4] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[2] +
    ked2_vsigma[0]*mgga_v4rhotau3[6] + mgga_v4rhosigmatau2[16]) + ked2_vsigma[0]*mgga_v4rho2tau2[7] +
    mgga_v4rho2sigmatau[16]) + ked2_vsigma[0]*mgga_v4rho2sigmatau[13] + mgga_v4rho2sigma2[14];
  v4rho2sigma2[15] = ked2_v2rho2[0]*mgga_v3sigma2tau[7] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigma2tau2[11] +
    2*ked2_vrho[0]*mgga_v4rhosigma2tau[19] + mgga_v4rho2sigma2[15];
  v4rho2sigma2[16] = ked2_v3rho2sigma[0]*mgga_v2sigmatau[3] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[5] +
    mgga_v3sigma2tau[9]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[7] + mgga_v4sigma2tau2[14]) +
    2*ked2_v2rhosigma[0]*mgga_v3rhosigmatau[9] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[5] +
    ked2_vsigma[0]*mgga_v4rhosigmatau2[14] + mgga_v4rhosigma2tau[21]) + ked2_vsigma[0]*mgga_v4rho2sigmatau[15] +
    mgga_v4rho2sigma2[16];
  v4rho2sigma2[17] = ked2_v4rho2sigma2[0]*mgga_vtau[1] + 2*ked2_v2rhosigma[0]*ked2_v2rhosigma[0]*mgga_v2tau2[2] +
    ked2_v2sigma2[0]*ked2_v2rho2[0]*mgga_v2tau2[2] + 2*ked2_vsigma[0]*ked2_v3rho2sigma[0]*mgga_v2tau2[2] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_v2rho2[0]*mgga_v3tau3[3] + 2*ked2_v3rho2sigma[0]*mgga_v2sigmatau[5] +
    2*ked2_vsigma[0]*ked2_v2rho2[0]*mgga_v3sigmatau2[8] + ked2_v2rho2[0]*mgga_v3sigma2tau[11] +
    ked2_vrho[0]*ked2_vrho[0]*(ked2_v2sigma2[0]*mgga_v3tau3[3] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[4] +
    2*ked2_vsigma[0]*mgga_v4sigmatau3[11] + mgga_v4sigma2tau2[17]) + 2*ked2_v3rhosigma2[0]*mgga_v2rhotau[3] +
    4*ked2_v2rhosigma[0]*(ked2_vrho[0]*mgga_v3sigmatau2[8] + ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[3] +
    mgga_v3rhotau2[5]) + mgga_v3rhosigmatau[11]) + 2*ked2_vrho[0]*(ked2_v3rhosigma2[0]*mgga_v2tau2[2] +
    ked2_v2sigma2[0]*mgga_v3rhotau2[5] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rhotau3[7] +
    2*ked2_vsigma[0]*mgga_v4rhosigmatau2[17] + mgga_v4rhosigma2tau[23]) + ked2_v2sigma2[0]*mgga_v3rho2tau[5] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rho2tau2[8] + 2*ked2_vsigma[0]*mgga_v4rho2sigmatau[17] +
    mgga_v4rho2sigma2[17];
  v4rho2sigmalapl[1] = ked1_v3rho2sigma[0]*mgga_v2lapltau[2] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3lapltau2[3] +
    mgga_v3sigmalapltau[2]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4lapltau3[4] + mgga_v4sigmalapltau2[3]) +
    2*ked1_v2rhosigma[0]*mgga_v3rholapltau[2] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3lapltau2[3] +
    ked1_vsigma[0]*mgga_v4rholapltau2[3] + mgga_v4rhosigmalapltau[2]) + ked1_vsigma[0]*mgga_v4rho2lapltau[2] +
    ked2_vlapl[0]*(ked1_v3rho2sigma[0]*mgga_v2tau2[1] + ked1_v2rho2[0]*(ked1_vsigma[0]*mgga_v3tau3[1] +
    mgga_v3sigmatau2[1]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) +
    2*ked1_v2rhosigma[0]*mgga_v3rhotau2[1] + 2*ked1_vrho[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[1] +
    ked1_vsigma[0]*mgga_v4rhotau3[1] + mgga_v4rhosigmatau2[1]) + ked1_vsigma[0]*mgga_v4rho2tau2[1] +
    mgga_v4rho2sigmatau[1]) + mgga_v4rho2sigmalapl[1];
  v4rho2sigmalapl[2] = ked1_v3rho2lapl[0]*mgga_v2sigmatau[2] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[3] +
    mgga_v3sigmalapltau[4]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[4] + mgga_v4sigmalapltau2[6]) +
    2*ked1_v2rholapl[0]*mgga_v3rhosigmatau[2] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[3] +
    ked1_vlapl[0]*mgga_v4rhosigmatau2[3] + mgga_v4rhosigmalapltau[4]) + ked1_vlapl[0]*mgga_v4rho2sigmatau[2] +
    mgga_v4rho2sigmalapl[2];
  v4rho2sigmalapl[3] = ked1_v2rho2[0]*mgga_v3sigmalapltau[6] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmalapltau2[9] +
    2*ked1_vrho[0]*mgga_v4rhosigmalapltau[6] + ked2_vlapl[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[5] + 2*ked1_vrho[0]*mgga_v4rhosigmatau2[4] + mgga_v4rho2sigmatau[3]) +
    mgga_v4rho2sigmalapl[3];
  v4rho2sigmalapl[4] = ked1_v3rho2lapl[0]*mgga_v2sigmatau[4] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[6] +
    mgga_v3sigmalapltau[8]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[8] + mgga_v4sigmalapltau2[12])
    + 2*ked1_v2rholapl[0]*mgga_v3rhosigmatau[4] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[6] +
    ked1_vlapl[0]*mgga_v4rhosigmatau2[6] + mgga_v4rhosigmalapltau[8]) +
    ked2_vsigma[0]*(ked1_v3rho2lapl[0]*mgga_v2tau2[1] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3tau3[1] +
    mgga_v3lapltau2[1]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) +
    2*ked1_v2rholapl[0]*mgga_v3rhotau2[1] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3tau3[1] +
    ked1_vlapl[0]*mgga_v4rhotau3[1] + mgga_v4rholapltau2[1]) + ked1_vlapl[0]*mgga_v4rho2tau2[1] +
    mgga_v4rho2lapltau[1]) + ked1_vlapl[0]*mgga_v4rho2sigmatau[4] + mgga_v4rho2sigmalapl[4];
  v4rho2sigmalapl[5] = ked2_vsigma[0]*ked1_v2rho2[0]*mgga_v3lapltau2[4] +
    ked2_vsigma[0]*ked1_vrho[0]*ked1_vrho[0]*mgga_v4lapltau3[5] + ked1_v2rho2[0]*mgga_v3sigmalapltau[10] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmalapltau2[15] + 2*ked2_vsigma[0]*ked1_vrho[0]*mgga_v4rholapltau2[4] +
    2*ked1_vrho[0]*mgga_v4rhosigmalapltau[10] + ked2_v2sigmalapl[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    ked2_vsigma[0]*mgga_v4rho2lapltau[3] + ked2_vlapl[0]*(ked1_v2rho2[0]*mgga_v3sigmatau2[7] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4sigmatau3[9] + 2*ked1_vrho[0]*mgga_v4rhosigmatau2[7] +
    ked2_vsigma[0]*(ked1_v2rho2[0]*mgga_v3tau3[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[2] +
    2*ked1_vrho[0]*mgga_v4rhotau3[2] + mgga_v4rho2tau2[2]) + mgga_v4rho2sigmatau[5]) + mgga_v4rho2sigmalapl[5];
  v4rho2sigmalapl[6] = ked1_v3rhosigmalapl[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    ked1_v2rhosigma[0]*(ked2_vrho[0]*mgga_v3lapltau2[1] + ked1_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[1] +
    mgga_v3rhotau2[3]) + mgga_v3rholapltau[4]) + ked1_v2rholapl[0]*(ked2_vrho[0]*mgga_v3sigmatau2[1] +
    ked1_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[3]) + mgga_v3rhosigmatau[6]) +
    ked1_vrho[0]*(ked2_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] + ked1_vsigma[0]*mgga_v4lapltau3[1] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) + mgga_v4sigmalapltau2[1]) +
    ked1_v2sigmalapl[0]*mgga_v3rhotau2[3] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4rhotau3[4] + mgga_v4rholapltau2[6]) +
    ked1_vlapl[0]*mgga_v4rhosigmatau2[9] + mgga_v4rhosigmalapltau[12]) +
    ked2_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3rhotau2[1] + ked1_vsigma[0]*mgga_v4rholapltau2[1] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4rhotau3[1] + mgga_v4rhosigmatau2[1]) + mgga_v4rhosigmalapltau[1]) +
    ked1_v2sigmalapl[0]*mgga_v3rho2tau[2] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4rho2tau2[3] + mgga_v4rho2lapltau[4]) +
    ked1_vlapl[0]*mgga_v4rho2sigmatau[6] + mgga_v4rho2sigmalapl[6];
  v4rho2sigmalapl[7] = ked1_v2rhosigma[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] + ked2_vrho[0]*mgga_v3lapltau2[4] +
    ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) + mgga_v3rholapltau[6]) +
    ked1_vrho[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[2] +
    mgga_v4sigmalapltau2[4]) + ked1_vsigma[0]*(ked2_v2rholapl[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4rhotau3[5] +
    mgga_v4rholapltau2[9]) + ked2_vlapl[0]*mgga_v4rhosigmatau2[10] + mgga_v4rhosigmalapltau[14]) +
    ked2_v2rholapl[0]*(ked1_vsigma[0]*mgga_v3rhotau2[1] + mgga_v3rhosigmatau[1]) +
    ked2_vrho[0]*(ked1_vsigma[0]*mgga_v4rholapltau2[4] + ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v4rhotau3[2] +
    mgga_v4rhosigmatau2[2]) + mgga_v4rhosigmalapltau[3]) + ked1_vsigma[0]*(ked2_vlapl[0]*mgga_v4rho2tau2[4] +
    mgga_v4rho2lapltau[6]) + ked2_vlapl[0]*mgga_v4rho2sigmatau[7] + mgga_v4rho2sigmalapl[7];
  v4rho2sigmalapl[8] = ked1_v2rholapl[0]*mgga_v3rhosigmatau[8] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[12] +
    mgga_v4rhosigmalapltau[16]) + ked2_vrho[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[4] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[5] + mgga_v4sigmalapltau2[7]) + ked1_vlapl[0]*mgga_v4rhosigmatau2[4] +
    mgga_v4rhosigmalapltau[5]) + ked1_vlapl[0]*mgga_v4rho2sigmatau[8] + mgga_v4rho2sigmalapl[8];
  v4rho2sigmalapl[9] = ked1_vrho[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[6] + mgga_v4sigmalapltau2[10]) + ked2_vlapl[0]*mgga_v4rhosigmatau2[13]
    + mgga_v4rhosigmalapltau[18]) + ked2_v2rholapl[0]*mgga_v3rhosigmatau[3] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[5] + mgga_v4rhosigmalapltau[7]) +
    ked2_vlapl[0]*mgga_v4rho2sigmatau[9] + mgga_v4rho2sigmalapl[9];
  v4rho2sigmalapl[10] = ked1_v2rholapl[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] + ked2_vrho[0]*mgga_v3sigmatau2[7] +
    ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) + mgga_v3rhosigmatau[10]) +
    ked1_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3lapltau2[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4lapltau3[2] +
    mgga_v4sigmalapltau2[13]) + ked2_vsigma[0]*mgga_v4rholapltau2[7] + ked1_vlapl[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[1]
    + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[2] + mgga_v4sigmatau3[9]) + ked2_vsigma[0]*mgga_v4rhotau3[5] +
    mgga_v4rhosigmatau2[15]) + mgga_v4rhosigmalapltau[20]) + ked2_v2rhosigma[0]*(ked1_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[1]) + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4rholapltau2[2] +
    ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4rhotau3[2] + mgga_v4rhosigmatau2[7]) + mgga_v4rhosigmalapltau[9]) +
    ked2_vsigma[0]*(ked1_vlapl[0]*mgga_v4rho2tau2[4] + mgga_v4rho2lapltau[5]) + ked1_vlapl[0]*mgga_v4rho2sigmatau[10] +
    mgga_v4rho2sigmalapl[10];
  v4rho2sigmalapl[11] = ked1_vrho[0]*(ked2_v3rhosigmalapl[0]*mgga_v2tau2[1] +
    ked2_v2rhosigma[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) +
    ked2_v2rholapl[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) +
    ked2_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] + ked2_vsigma[0]*mgga_v4lapltau3[6] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) + mgga_v4sigmalapltau2[16]) +
    ked2_v2sigmalapl[0]*mgga_v3rhotau2[4] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4rhotau3[6] + mgga_v4rholapltau2[10]) +
    ked2_vlapl[0]*mgga_v4rhosigmatau2[16] + mgga_v4rhosigmalapltau[22]) + ked2_v3rhosigmalapl[0]*mgga_v2rhotau[1] +
    ked2_v2rhosigma[0]*(ked2_vlapl[0]*mgga_v3rhotau2[2] + mgga_v3rholapltau[3]) +
    ked2_v2rholapl[0]*(ked2_vsigma[0]*mgga_v3rhotau2[2] + mgga_v3rhosigmatau[5]) +
    ked2_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3rhotau2[2] + ked2_vsigma[0]*mgga_v4rholapltau2[5] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4rhotau3[3] + mgga_v4rhosigmatau2[8]) + mgga_v4rhosigmalapltau[11]) +
    ked2_v2sigmalapl[0]*mgga_v3rho2tau[3] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4rho2tau2[5] + mgga_v4rho2lapltau[7]) +
    ked2_vlapl[0]*mgga_v4rho2sigmatau[11] + mgga_v4rho2sigmalapl[11];
  v4rho2sigmalapl[12] = ked1_vsigma[0]*ked2_v2rho2[0]*mgga_v3lapltau2[1] +
    ked1_vsigma[0]*ked2_vrho[0]*ked2_vrho[0]*mgga_v4lapltau3[2] + ked2_v2rho2[0]*mgga_v3sigmalapltau[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmalapltau2[2] + 2*ked1_vsigma[0]*ked2_vrho[0]*mgga_v4rholapltau2[7] +
    2*ked2_vrho[0]*mgga_v4rhosigmalapltau[13] + ked1_v2sigmalapl[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] + 2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) +
    ked1_vsigma[0]*mgga_v4rho2lapltau[8] + ked1_vlapl[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[2] + 2*ked2_vrho[0]*mgga_v4rhosigmatau2[10] +
    ked1_vsigma[0]*(ked2_v2rho2[0]*mgga_v3tau3[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[2] +
    2*ked2_vrho[0]*mgga_v4rhotau3[5] + mgga_v4rho2tau2[6]) + mgga_v4rho2sigmatau[12]) + mgga_v4rho2sigmalapl[12];
  v4rho2sigmalapl[13] = ked2_v3rho2lapl[0]*mgga_v2sigmatau[1] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[2] +
    mgga_v3sigmalapltau[3]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[3] + mgga_v4sigmalapltau2[5]) +
    2*ked2_v2rholapl[0]*mgga_v3rhosigmatau[7] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[2] +
    ked2_vlapl[0]*mgga_v4rhosigmatau2[11] + mgga_v4rhosigmalapltau[15]) +
    ked1_vsigma[0]*(ked2_v3rho2lapl[0]*mgga_v2tau2[1] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3tau3[2] +
    mgga_v3lapltau2[4]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[3] + mgga_v4lapltau3[6]) +
    2*ked2_v2rholapl[0]*mgga_v3rhotau2[4] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3tau3[2] +
    ked2_vlapl[0]*mgga_v4rhotau3[6] + mgga_v4rholapltau2[10]) + ked2_vlapl[0]*mgga_v4rho2tau2[7] +
    mgga_v4rho2lapltau[10]) + ked2_vlapl[0]*mgga_v4rho2sigmatau[13] + mgga_v4rho2sigmalapl[13];
  v4rho2sigmalapl[14] = ked2_v2rho2[0]*mgga_v3sigmalapltau[5] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmalapltau2[8] +
    2*ked2_vrho[0]*mgga_v4rhosigmalapltau[17] + ked1_vlapl[0]*(ked2_v2rho2[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4sigmatau3[6] + 2*ked2_vrho[0]*mgga_v4rhosigmatau2[13] + mgga_v4rho2sigmatau[14]) +
    mgga_v4rho2sigmalapl[14];
  v4rho2sigmalapl[15] = ked2_v3rho2lapl[0]*mgga_v2sigmatau[3] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[5] +
    mgga_v3sigmalapltau[7]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[7] + mgga_v4sigmalapltau2[11])
    + 2*ked2_v2rholapl[0]*mgga_v3rhosigmatau[9] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[5] +
    ked2_vlapl[0]*mgga_v4rhosigmatau2[14] + mgga_v4rhosigmalapltau[19]) + ked2_vlapl[0]*mgga_v4rho2sigmatau[15] +
    mgga_v4rho2sigmalapl[15];
  v4rho2sigmalapl[16] = ked2_v3rho2sigma[0]*mgga_v2lapltau[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3lapltau2[2] +
    mgga_v3sigmalapltau[9]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4lapltau3[3] + mgga_v4sigmalapltau2[14])
    + 2*ked2_v2rhosigma[0]*mgga_v3rholapltau[5] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3lapltau2[2] +
    ked2_vsigma[0]*mgga_v4rholapltau2[8] + mgga_v4rhosigmalapltau[21]) + ked2_vsigma[0]*mgga_v4rho2lapltau[9] +
    ked1_vlapl[0]*(ked2_v3rho2sigma[0]*mgga_v2tau2[1] + ked2_v2rho2[0]*(ked2_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[7]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) +
    2*ked2_v2rhosigma[0]*mgga_v3rhotau2[4] + 2*ked2_vrho[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[2] +
    ked2_vsigma[0]*mgga_v4rhotau3[6] + mgga_v4rhosigmatau2[16]) + ked2_vsigma[0]*mgga_v4rho2tau2[7] +
    mgga_v4rho2sigmatau[16]) + mgga_v4rho2sigmalapl[16];
  v4rho2sigmalapl[17] = ked2_v4rho2sigmalapl[0]*mgga_vtau[1] + ked2_v2sigmalapl[0]*ked2_v2rho2[0]*mgga_v2tau2[2] +
    ked2_vsigma[0]*ked2_v3rho2lapl[0]*mgga_v2tau2[2] + ked2_vlapl[0]*ked2_v3rho2sigma[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*ked2_vsigma[0]*ked2_v2rho2[0]*mgga_v3tau3[3] + ked2_v3rho2sigma[0]*mgga_v2lapltau[3] +
    ked2_vsigma[0]*ked2_v2rho2[0]*mgga_v3lapltau2[5] + ked2_v3rho2lapl[0]*mgga_v2sigmatau[5] +
    ked2_vlapl[0]*ked2_v2rho2[0]*mgga_v3sigmatau2[8] + ked2_v2rho2[0]*mgga_v3sigmalapltau[11] +
    ked2_vrho[0]*ked2_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[3] + ked2_vsigma[0]*mgga_v4lapltau3[7] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[4] + mgga_v4sigmatau3[11]) + mgga_v4sigmalapltau2[17]) +
    2*ked2_v3rhosigmalapl[0]*mgga_v2rhotau[3] + 2*ked2_vlapl[0]*ked2_v2rhosigma[0]*mgga_v3rhotau2[5] +
    2*ked2_v2rhosigma[0]*mgga_v3rholapltau[7] + 2*ked2_v2rholapl[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[2] +
    ked2_vrho[0]*mgga_v3sigmatau2[8] + ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[3] + mgga_v3rhotau2[5]) +
    mgga_v3rhosigmatau[11]) + 2*ked2_vrho[0]*(ked2_v3rhosigmalapl[0]*mgga_v2tau2[2] +
    ked2_v2rhosigma[0]*mgga_v3lapltau2[5] + ked2_v2sigmalapl[0]*mgga_v3rhotau2[5] +
    ked2_vsigma[0]*mgga_v4rholapltau2[11] + ked2_vlapl[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[3] +
    ked2_vsigma[0]*mgga_v4rhotau3[7] + mgga_v4rhosigmatau2[17]) + mgga_v4rhosigmalapltau[23]) +
    ked2_v2sigmalapl[0]*mgga_v3rho2tau[5] + ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4rho2tau2[8] +
    ked2_vsigma[0]*mgga_v4rho2lapltau[11] + ked2_vlapl[0]*mgga_v4rho2sigmatau[17] + mgga_v4rho2sigmalapl[17];
  v4rho2sigmatau[1] = 0;
  v4rho2sigmatau[2] = 0;
  v4rho2sigmatau[3] = 0;
  v4rho2sigmatau[4] = 0;
  v4rho2sigmatau[5] = 0;
  v4rho2sigmatau[6] = 0;
  v4rho2sigmatau[7] = 0;
  v4rho2sigmatau[8] = 0;
  v4rho2sigmatau[9] = 0;
  v4rho2sigmatau[10] = 0;
  v4rho2sigmatau[11] = 0;
  v4rho2sigmatau[12] = 0;
  v4rho2sigmatau[13] = 0;
  v4rho2sigmatau[14] = 0;
  v4rho2sigmatau[15] = 0;
  v4rho2sigmatau[16] = 0;
  v4rho2sigmatau[17] = 0;
  v4rho2lapl2[1] = ked1_v3rho2lapl[0]*mgga_v2lapltau[2] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3lapltau2[3] +
    mgga_v3lapl2tau[2]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4lapltau3[4] + mgga_v4lapl2tau2[3]) +
    2*ked1_v2rholapl[0]*mgga_v3rholapltau[2] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3lapltau2[3] +
    ked1_vlapl[0]*mgga_v4rholapltau2[3] + mgga_v4rholapl2tau[2]) + ked1_vlapl[0]*mgga_v4rho2lapltau[2] +
    ked2_vlapl[0]*(ked1_v3rho2lapl[0]*mgga_v2tau2[1] + ked1_v2rho2[0]*(ked1_vlapl[0]*mgga_v3tau3[1] +
    mgga_v3lapltau2[1]) + ked1_vrho[0]*ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) +
    2*ked1_v2rholapl[0]*mgga_v3rhotau2[1] + 2*ked1_vrho[0]*(ked1_v2rholapl[0]*mgga_v3tau3[1] +
    ked1_vlapl[0]*mgga_v4rhotau3[1] + mgga_v4rholapltau2[1]) + ked1_vlapl[0]*mgga_v4rho2tau2[1] +
    mgga_v4rho2lapltau[1]) + mgga_v4rho2lapl2[1];
  v4rho2lapl2[2] = ked1_v2rho2[0]*mgga_v3lapl2tau[4] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4lapl2tau2[6] +
    2*ked1_vrho[0]*mgga_v4rholapl2tau[4] + ked2_v2lapl2[0]*(ked1_v2rho2[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v3tau3[1] + 2*ked1_vrho[0]*mgga_v3rhotau2[1] + mgga_v3rho2tau[1]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked1_v2rho2[0]*mgga_v3tau3[2] + ked1_vrho[0]*ked1_vrho[0]*mgga_v4tau4[2] +
    2*ked1_vrho[0]*mgga_v4rhotau3[2] + mgga_v4rho2tau2[2]) + 2*ked2_vlapl[0]*(ked1_v2rho2[0]*mgga_v3lapltau2[4] +
    ked1_vrho[0]*ked1_vrho[0]*mgga_v4lapltau3[5] + 2*ked1_vrho[0]*mgga_v4rholapltau2[4] + mgga_v4rho2lapltau[3]) +
    mgga_v4rho2lapl2[2];
  v4rho2lapl2[3] = ked1_v3rholapl2[0]*mgga_v2rhotau[2] + 2*ked1_vlapl[0]*ked1_v2rholapl[0]*mgga_v3rhotau2[3] +
    ked1_vlapl[0]*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4rhotau3[4] + 2*ked1_v2rholapl[0]*mgga_v3rholapltau[4] +
    2*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4rholapltau2[6] + ked1_vrho[0]*mgga_v4rholapl2tau[6] +
    ked2_vrho[0]*(ked1_v3rholapl2[0]*mgga_v2tau2[1] + 2*ked1_v2rholapl[0]*mgga_v3lapltau2[1] +
    ked1_vrho[0]*mgga_v4lapl2tau2[1] + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[1]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[1]) +
    2*ked1_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3tau3[1] + ked1_vrho[0]*mgga_v4lapltau3[1] + mgga_v4rholapltau2[1]) +
    mgga_v4rholapl2tau[1]) + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3rhotau2[3] + mgga_v3rho2tau[2]) +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rho2tau2[3] + 2*ked1_vlapl[0]*mgga_v4rho2lapltau[4] + mgga_v4rho2lapl2[3];
  v4rho2lapl2[4] = ked1_v2rholapl[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] + ked2_vrho[0]*mgga_v3lapltau2[4] +
    ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) + mgga_v3rholapltau[6]) +
    ked1_vrho[0]*(ked2_v2rholapl[0]*mgga_v3lapltau2[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4lapltau3[2] +
    mgga_v4lapl2tau2[4]) + ked1_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3tau3[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[2]
    + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4rhotau3[5] + mgga_v4rholapltau2[9]) +
    ked2_vlapl[0]*mgga_v4rholapltau2[7] + mgga_v4rholapl2tau[8]) + ked2_v2rholapl[0]*(ked1_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[1]) + ked2_vrho[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhotau3[2] + mgga_v4rholapltau2[4]) +
    ked2_vlapl[0]*mgga_v4rholapltau2[2] + mgga_v4rholapl2tau[3]) + ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rho2tau2[4] +
    mgga_v4rho2lapltau[6]) + ked2_vlapl[0]*mgga_v4rho2lapltau[5] + mgga_v4rho2lapl2[4];
  v4rho2lapl2[5] = ked1_vrho[0]*(ked2_v3rholapl2[0]*mgga_v2tau2[1] + 2*ked2_v2rholapl[0]*mgga_v3lapltau2[4] +
    ked2_vrho[0]*mgga_v4lapl2tau2[7] + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[6]) +
    2*ked2_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3tau3[2] + ked2_vrho[0]*mgga_v4lapltau3[6] + mgga_v4rholapltau2[10]) +
    mgga_v4rholapl2tau[10]) + ked2_v3rholapl2[0]*mgga_v2rhotau[1] + 2*ked2_vlapl[0]*ked2_v2rholapl[0]*mgga_v3rhotau2[2]
    + ked2_vlapl[0]*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4rhotau3[3] + 2*ked2_v2rholapl[0]*mgga_v3rholapltau[3] +
    2*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4rholapltau2[5] + ked2_vrho[0]*mgga_v4rholapl2tau[5] +
    ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3rhotau2[2] + mgga_v3rho2tau[3]) +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rho2tau2[5] + 2*ked2_vlapl[0]*mgga_v4rho2lapltau[7] + mgga_v4rho2lapl2[5];
  v4rho2lapl2[6] = ked2_v2rho2[0]*mgga_v3lapl2tau[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4lapl2tau2[2] +
    2*ked2_vrho[0]*mgga_v4rholapl2tau[7] + ked1_v2lapl2[0]*(ked2_v2rho2[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v3tau3[2] + 2*ked2_vrho[0]*mgga_v3rhotau2[4] + mgga_v3rho2tau[4]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked2_v2rho2[0]*mgga_v3tau3[1] + ked2_vrho[0]*ked2_vrho[0]*mgga_v4tau4[2] +
    2*ked2_vrho[0]*mgga_v4rhotau3[5] + mgga_v4rho2tau2[6]) + 2*ked1_vlapl[0]*(ked2_v2rho2[0]*mgga_v3lapltau2[1] +
    ked2_vrho[0]*ked2_vrho[0]*mgga_v4lapltau3[2] + 2*ked2_vrho[0]*mgga_v4rholapltau2[7] + mgga_v4rho2lapltau[8]) +
    mgga_v4rho2lapl2[6];
  v4rho2lapl2[7] = ked2_v3rho2lapl[0]*mgga_v2lapltau[1] + ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3lapltau2[2] +
    mgga_v3lapl2tau[3]) + ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4lapltau3[3] + mgga_v4lapl2tau2[5]) +
    2*ked2_v2rholapl[0]*mgga_v3rholapltau[5] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3lapltau2[2] +
    ked2_vlapl[0]*mgga_v4rholapltau2[8] + mgga_v4rholapl2tau[9]) + ked1_vlapl[0]*(ked2_v3rho2lapl[0]*mgga_v2tau2[1] +
    ked2_v2rho2[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) +
    ked2_vrho[0]*ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[3] + mgga_v4lapltau3[6]) +
    2*ked2_v2rholapl[0]*mgga_v3rhotau2[4] + 2*ked2_vrho[0]*(ked2_v2rholapl[0]*mgga_v3tau3[2] +
    ked2_vlapl[0]*mgga_v4rhotau3[6] + mgga_v4rholapltau2[10]) + ked2_vlapl[0]*mgga_v4rho2tau2[7] +
    mgga_v4rho2lapltau[10]) + ked2_vlapl[0]*mgga_v4rho2lapltau[9] + mgga_v4rho2lapl2[7];
  v4rho2lapl2[8] = ked2_v4rho2lapl2[0]*mgga_vtau[1] + 2*ked2_v2rholapl[0]*ked2_v2rholapl[0]*mgga_v2tau2[2] +
    ked2_v2lapl2[0]*ked2_v2rho2[0]*mgga_v2tau2[2] + 2*ked2_vlapl[0]*ked2_v3rho2lapl[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_v2rho2[0]*mgga_v3tau3[3] + 2*ked2_v3rho2lapl[0]*mgga_v2lapltau[3] +
    2*ked2_vlapl[0]*ked2_v2rho2[0]*mgga_v3lapltau2[5] + ked2_v2rho2[0]*mgga_v3lapl2tau[5] +
    ked2_vrho[0]*ked2_vrho[0]*(ked2_v2lapl2[0]*mgga_v3tau3[3] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4tau4[4] +
    2*ked2_vlapl[0]*mgga_v4lapltau3[7] + mgga_v4lapl2tau2[8]) + 2*ked2_v3rholapl2[0]*mgga_v2rhotau[3] +
    4*ked2_v2rholapl[0]*(ked2_vrho[0]*mgga_v3lapltau2[5] + ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[3] +
    mgga_v3rhotau2[5]) + mgga_v3rholapltau[7]) + 2*ked2_vrho[0]*(ked2_v3rholapl2[0]*mgga_v2tau2[2] +
    ked2_v2lapl2[0]*mgga_v3rhotau2[5] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rhotau3[7] +
    2*ked2_vlapl[0]*mgga_v4rholapltau2[11] + mgga_v4rholapl2tau[11]) + ked2_v2lapl2[0]*mgga_v3rho2tau[5] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rho2tau2[8] + 2*ked2_vlapl[0]*mgga_v4rho2lapltau[11] + mgga_v4rho2lapl2[8];
  v4rho2lapltau[1] = 0;
  v4rho2lapltau[2] = 0;
  v4rho2lapltau[3] = 0;
  v4rho2lapltau[4] = 0;
  v4rho2lapltau[5] = 0;
  v4rho2lapltau[6] = 0;
  v4rho2lapltau[7] = 0;
  v4rho2lapltau[8] = 0;
  v4rho2lapltau[9] = 0;
  v4rho2lapltau[10] = 0;
  v4rho2lapltau[11] = 0;
  v4rho2tau2[1] = 0;
  v4rho2tau2[2] = 0;
  v4rho2tau2[3] = 0;
  v4rho2tau2[4] = 0;
  v4rho2tau2[5] = 0;
  v4rho2tau2[6] = 0;
  v4rho2tau2[7] = 0;
  v4rho2tau2[8] = 0;
  v4rhosigma3[1] = ked1_v3rhosigma2[0]*mgga_v2sigmatau[2] + 2*ked1_v2rhosigma[0]*mgga_v3sigma2tau[2] +
    ked1_vrho[0]*mgga_v4sigma3tau[2] + ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3sigmatau2[3] + mgga_v3rhosigmatau[2]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vrho[0]*mgga_v4sigmatau3[4] + mgga_v4rhosigmatau2[3]) +
    2*ked1_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[3] + ked1_vrho[0]*mgga_v4sigma2tau2[3] +
    mgga_v4rhosigma2tau[2]) + mgga_v4rhosigma3[1];
  v4rhosigma3[2] = ked1_v3rhosigma2[0]*mgga_v2sigmatau[4] + 2*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3sigmatau2[6] +
    ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigmatau3[8] + 2*ked1_v2rhosigma[0]*mgga_v3sigma2tau[4] +
    2*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigma2tau2[6] + ked1_vrho[0]*mgga_v4sigma3tau[4] +
    ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3sigmatau2[6] + mgga_v3rhosigmatau[4]) +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rhosigmatau2[6] + 2*ked1_vsigma[0]*mgga_v4rhosigma2tau[4] +
    ked2_vsigma[0]*(ked1_v3rhosigma2[0]*mgga_v2tau2[1] + 2*ked1_v2rhosigma[0]*mgga_v3sigmatau2[1] +
    ked1_vrho[0]*mgga_v4sigma2tau2[1] + ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[1]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[1]) +
    2*ked1_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[1] + ked1_vrho[0]*mgga_v4sigmatau3[1] + mgga_v4rhosigmatau2[1]) +
    mgga_v4rhosigma2tau[1]) + mgga_v4rhosigma3[2];
  v4rhosigma3[3] = ked1_v2rhosigma[0]*mgga_v3sigma2tau[6] + ked1_vrho[0]*mgga_v4sigma3tau[6] +
    ked1_vsigma[0]*(ked1_vrho[0]*mgga_v4sigma2tau2[9] + mgga_v4rhosigma2tau[6]) + mgga_v4rhosigma3[3];
  v4rhosigma3[4] = ked1_v2rhosigma[0]*mgga_v3sigma2tau[8] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigma2tau2[12] +
    mgga_v4sigma3tau[8]) + ked1_vsigma[0]*mgga_v4rhosigma2tau[8] +
    ked2_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[4] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[5] +
    mgga_v4sigma2tau2[4]) + ked1_vsigma[0]*mgga_v4rhosigmatau2[4] + mgga_v4rhosigma2tau[3]) + mgga_v4rhosigma3[4];
  v4rhosigma3[5] = ked1_v2rhosigma[0]*mgga_v3sigma2tau[10] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigma2tau2[15] +
    mgga_v4sigma3tau[10]) + ked1_vsigma[0]*mgga_v4rhosigma2tau[10] +
    ked2_v2sigma2[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[1] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3tau3[1] +
    mgga_v3sigmatau2[1]) + ked1_vsigma[0]*mgga_v3rhotau2[1] + mgga_v3rhosigmatau[1]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[2] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4tau4[2] +
    mgga_v4sigmatau3[2]) + ked1_vsigma[0]*mgga_v4rhotau3[2] + mgga_v4rhosigmatau2[2]) +
    2*ked2_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[7] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[9] +
    mgga_v4sigma2tau2[7]) + ked1_vsigma[0]*mgga_v4rhosigmatau2[7] + mgga_v4rhosigma2tau[5]) + mgga_v4rhosigma3[5];
  v4rhosigma3[6] = ked1_vrho[0]*mgga_v4sigma3tau[12] + mgga_v4rhosigma3[6];
  v4rhosigma3[7] = ked1_vrho[0]*mgga_v4sigma3tau[14] + ked2_vsigma[0]*(ked1_vrho[0]*mgga_v4sigma2tau2[10] +
    mgga_v4rhosigma2tau[7]) + mgga_v4rhosigma3[7];
  v4rhosigma3[8] = ked1_vrho[0]*mgga_v4sigma3tau[16] + ked2_v2sigma2[0]*(ked1_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[3]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v4sigmatau3[6] + mgga_v4rhosigmatau2[5])
    + 2*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v4sigma2tau2[13] + mgga_v4rhosigma2tau[9]) + mgga_v4rhosigma3[8];
  v4rhosigma3[9] = ked1_vrho[0]*(3*ked2_v2sigma2[0]*mgga_v3sigmatau2[7] + mgga_v4sigma3tau[18]) +
    ked2_v3sigma3[0]*(ked1_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[1]) +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[3]) +
    3*ked2_v2sigma2[0]*mgga_v3rhosigmatau[5] + 3*ked2_vsigma[0]*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v4sigmatau3[10] +
    mgga_v4rhosigmatau2[8]) + 3*ked2_vsigma[0]*(ked1_vrho[0]*(ked2_v2sigma2[0]*mgga_v3tau3[2] + mgga_v4sigma2tau2[16])
    + ked2_v2sigma2[0]*mgga_v3rhotau2[2] + mgga_v4rhosigma2tau[11]) + mgga_v4rhosigma3[9];
  v4rhosigma3[10] = ked2_vrho[0]*(3*ked1_v2sigma2[0]*mgga_v3sigmatau2[1] + mgga_v4sigma3tau[1]) +
    ked1_v3sigma3[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[4]) +
    3*ked1_v2sigma2[0]*mgga_v3rhosigmatau[6] + 3*ked1_vsigma[0]*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v4sigmatau3[1] +
    mgga_v4rhosigmatau2[9]) + 3*ked1_vsigma[0]*(ked2_vrho[0]*(ked1_v2sigma2[0]*mgga_v3tau3[1] + mgga_v4sigma2tau2[1]) +
    ked1_v2sigma2[0]*mgga_v3rhotau2[3] + mgga_v4rhosigma2tau[12]) + mgga_v4rhosigma3[10];
  v4rhosigma3[11] = ked2_vrho[0]*mgga_v4sigma3tau[3] + ked1_v2sigma2[0]*(ked2_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[8]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v4sigmatau3[5] + mgga_v4rhosigmatau2[12])
    + 2*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v4sigma2tau2[4] + mgga_v4rhosigma2tau[14]) + mgga_v4rhosigma3[11];
  v4rhosigma3[12] = ked2_v2rhosigma[0]*mgga_v3sigma2tau[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigma2tau2[2] +
    mgga_v4sigma3tau[5]) + ked1_v2sigma2[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + ked2_vsigma[0]*mgga_v3rhotau2[4] +
    mgga_v3rhosigmatau[10]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[2] + mgga_v4sigmatau3[9]) + ked2_vsigma[0]*mgga_v4rhotau3[5] +
    mgga_v4rhosigmatau2[15]) + 2*ked1_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[1] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[2] + mgga_v4sigma2tau2[7]) + ked2_vsigma[0]*mgga_v4rhosigmatau2[10] +
    mgga_v4rhosigma2tau[16]) + ked2_vsigma[0]*mgga_v4rhosigma2tau[13] + mgga_v4rhosigma3[12];
  v4rhosigma3[13] = ked2_vrho[0]*mgga_v4sigma3tau[7] + ked1_vsigma[0]*(ked2_vrho[0]*mgga_v4sigma2tau2[10] +
    mgga_v4rhosigma2tau[18]) + mgga_v4rhosigma3[13];
  v4rhosigma3[14] = ked2_v2rhosigma[0]*mgga_v3sigma2tau[3] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigma2tau2[5] +
    mgga_v4sigma3tau[9]) + ked1_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[6] + mgga_v4sigma2tau2[13]) + ked2_vsigma[0]*mgga_v4rhosigmatau2[13]
    + mgga_v4rhosigma2tau[20]) + ked2_vsigma[0]*mgga_v4rhosigma2tau[15] + mgga_v4rhosigma3[14];
  v4rhosigma3[15] = ked2_v3rhosigma2[0]*mgga_v2sigmatau[1] + 2*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3sigmatau2[2] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigmatau3[3] + 2*ked2_v2rhosigma[0]*mgga_v3sigma2tau[5] +
    2*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigma2tau2[8] + ked2_vrho[0]*mgga_v4sigma3tau[11] +
    ked1_vsigma[0]*(ked2_v3rhosigma2[0]*mgga_v2tau2[1] + 2*ked2_v2rhosigma[0]*mgga_v3sigmatau2[7] +
    ked2_vrho[0]*mgga_v4sigma2tau2[16] + ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[6]) +
    2*ked2_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[2] + ked2_vrho[0]*mgga_v4sigmatau3[10] + mgga_v4rhosigmatau2[16])
    + mgga_v4rhosigma2tau[22]) + ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3sigmatau2[2] + mgga_v3rhosigmatau[7]) +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rhosigmatau2[11] + 2*ked2_vsigma[0]*mgga_v4rhosigma2tau[17] +
    mgga_v4rhosigma3[15];
  v4rhosigma3[16] = ked2_vrho[0]*mgga_v4sigma3tau[13] + mgga_v4rhosigma3[16];
  v4rhosigma3[17] = ked2_v2rhosigma[0]*mgga_v3sigma2tau[7] + ked2_vrho[0]*mgga_v4sigma3tau[15] +
    ked2_vsigma[0]*(ked2_vrho[0]*mgga_v4sigma2tau2[11] + mgga_v4rhosigma2tau[19]) + mgga_v4rhosigma3[17];
  v4rhosigma3[18] = ked2_v3rhosigma2[0]*mgga_v2sigmatau[3] + 2*ked2_v2rhosigma[0]*mgga_v3sigma2tau[9] +
    ked2_vrho[0]*mgga_v4sigma3tau[17] + ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3sigmatau2[5] + mgga_v3rhosigmatau[9]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vrho[0]*mgga_v4sigmatau3[7] + mgga_v4rhosigmatau2[14]) +
    2*ked2_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[5] + ked2_vrho[0]*mgga_v4sigma2tau2[14] +
    mgga_v4rhosigma2tau[21]) + mgga_v4rhosigma3[18];
  v4rhosigma3[19] = ked2_v4rhosigma3[0]*mgga_vtau[1] + 3*ked2_vsigma[0]*ked2_v3rhosigma2[0]*mgga_v2tau2[2] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3tau3[3] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4tau4[4] + 3*ked2_v3rhosigma2[0]*mgga_v2sigmatau[5]
    + 6*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3sigmatau2[8] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigmatau3[11] + 3*ked2_v2rhosigma[0]*mgga_v3sigma2tau[11] +
    3*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigma2tau2[17] + ked2_vrho[0]*mgga_v4sigma3tau[19] +
    ked2_v3sigma3[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rhotau3[7] +
    3*ked2_v2sigma2[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[2] + ked2_vrho[0]*mgga_v3sigmatau2[8] +
    ked2_vsigma[0]*(ked2_vrho[0]*mgga_v3tau3[3] + mgga_v3rhotau2[5]) + mgga_v3rhosigmatau[11]) +
    3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rhosigmatau2[17] + 3*ked2_vsigma[0]*mgga_v4rhosigma2tau[23] +
    mgga_v4rhosigma3[19];
  v4rhosigma2lapl[1] = ked1_v3rhosigma2[0]*mgga_v2lapltau[2] + 2*ked1_vsigma[0]*ked1_v2rhosigma[0]*mgga_v3lapltau2[3] +
    ked1_vsigma[0]*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4lapltau3[4] + 2*ked1_v2rhosigma[0]*mgga_v3sigmalapltau[2] +
    2*ked1_vsigma[0]*ked1_vrho[0]*mgga_v4sigmalapltau2[3] + ked1_vrho[0]*mgga_v4sigma2lapltau[2] +
    ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3lapltau2[3] + mgga_v3rholapltau[2]) +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4rholapltau2[3] + 2*ked1_vsigma[0]*mgga_v4rhosigmalapltau[2] +
    ked2_vlapl[0]*(ked1_v3rhosigma2[0]*mgga_v2tau2[1] + 2*ked1_v2rhosigma[0]*mgga_v3sigmatau2[1] +
    ked1_vrho[0]*mgga_v4sigma2tau2[1] + ked1_v2sigma2[0]*(ked1_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[1]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[1]) +
    2*ked1_vsigma[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[1] + ked1_vrho[0]*mgga_v4sigmatau3[1] + mgga_v4rhosigmatau2[1]) +
    mgga_v4rhosigma2tau[1]) + mgga_v4rhosigma2lapl[1];
  v4rhosigma2lapl[2] = ked1_v3rhosigmalapl[0]*mgga_v2sigmatau[2] +
    ked1_v2rhosigma[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[3] + mgga_v3sigmalapltau[4]) +
    ked1_v2rholapl[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[3] + mgga_v3sigma2tau[2]) +
    ked1_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[3] + ked1_vsigma[0]*mgga_v4sigmalapltau2[6] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[4] + mgga_v4sigma2tau2[3]) + mgga_v4sigma2lapltau[4]) +
    ked1_v2sigmalapl[0]*mgga_v3rhosigmatau[2] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[3] +
    mgga_v4rhosigmalapltau[4]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[2] + mgga_v4rhosigma2lapl[2];
  v4rhosigma2lapl[3] = ked1_v2rhosigma[0]*mgga_v3sigmalapltau[6] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmalapltau2[9]
    + mgga_v4sigma2lapltau[6]) + ked1_vsigma[0]*mgga_v4rhosigmalapltau[6] +
    ked2_vlapl[0]*(ked1_v2rhosigma[0]*mgga_v3sigmatau2[4] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[5] +
    mgga_v4sigma2tau2[4]) + ked1_vsigma[0]*mgga_v4rhosigmatau2[4] + mgga_v4rhosigma2tau[3]) + mgga_v4rhosigma2lapl[3];
  v4rhosigma2lapl[4] = ked1_v3rhosigmalapl[0]*(ked2_vsigma[0]*mgga_v2tau2[1] + mgga_v2sigmatau[4]) +
    ked1_v2rhosigma[0]*(ked2_vsigma[0]*mgga_v3lapltau2[1] + ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[1] +
    mgga_v3sigmatau2[6]) + mgga_v3sigmalapltau[8]) + ked1_v2rholapl[0]*(ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v3tau3[1] +
    mgga_v3sigmatau2[6]) + ked2_vsigma[0]*mgga_v3sigmatau2[1] + mgga_v3sigma2tau[4]) +
    ked1_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[6] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[8] +
    mgga_v4sigmalapltau2[12]) + ked2_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) + ked1_vlapl[0]*mgga_v4sigmatau3[1] +
    mgga_v4sigmalapltau2[1]) + ked1_vlapl[0]*mgga_v4sigma2tau2[6] + mgga_v4sigma2lapltau[8]) +
    ked1_v2sigmalapl[0]*(ked2_vsigma[0]*mgga_v3rhotau2[1] + mgga_v3rhosigmatau[4]) +
    ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v4rholapltau2[1] + ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4rhotau3[1] +
    mgga_v4rhosigmatau2[6]) + mgga_v4rhosigmalapltau[8]) + ked2_vsigma[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[1] +
    mgga_v4rhosigmalapltau[1]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[4] + mgga_v4rhosigma2lapl[4];
  v4rhosigma2lapl[5] = ked1_v2rhosigma[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] + ked2_vsigma[0]*mgga_v3lapltau2[4] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + mgga_v3sigmalapltau[10]) +
    ked1_vrho[0]*(ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[2] +
    mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4sigmatau3[9] + mgga_v4sigmalapltau2[15]) +
    ked2_v2sigmalapl[0]*mgga_v3sigmatau2[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[2] +
    mgga_v4sigmalapltau2[4]) + ked2_vlapl[0]*mgga_v4sigma2tau2[7] + mgga_v4sigma2lapltau[10]) +
    ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3rhotau2[1] + ked2_vsigma[0]*mgga_v4rholapltau2[4] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4rhotau3[2] + mgga_v4rhosigmatau2[7]) + mgga_v4rhosigmalapltau[10]) +
    ked2_v2sigmalapl[0]*mgga_v3rhosigmatau[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[2] +
    mgga_v4rhosigmalapltau[3]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[5] + mgga_v4rhosigma2lapl[5];
  v4rhosigma2lapl[6] = ked1_v2rholapl[0]*mgga_v3sigma2tau[6] + ked1_vrho[0]*mgga_v4sigma2lapltau[12] +
    ked1_vlapl[0]*(ked1_vrho[0]*mgga_v4sigma2tau2[9] + mgga_v4rhosigma2tau[6]) + mgga_v4rhosigma2lapl[6];
  v4rhosigma2lapl[7] = ked1_vrho[0]*mgga_v4sigma2lapltau[14] + ked2_vlapl[0]*(ked1_vrho[0]*mgga_v4sigma2tau2[10] +
    mgga_v4rhosigma2tau[7]) + mgga_v4rhosigma2lapl[7];
  v4rhosigma2lapl[8] = ked1_v2rholapl[0]*mgga_v3sigma2tau[8] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigma2tau2[12] +
    mgga_v4sigma2lapltau[16]) + ked2_vsigma[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[4] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[5] + mgga_v4sigmalapltau2[7]) + ked1_vlapl[0]*mgga_v4rhosigmatau2[4] +
    mgga_v4rhosigmalapltau[5]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[8] + mgga_v4rhosigma2lapl[8];
  v4rhosigma2lapl[9] = ked1_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[4] + ked2_vsigma[0]*mgga_v4sigmalapltau2[10]
    + ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[6] + mgga_v4sigma2tau2[13]) + mgga_v4sigma2lapltau[18]) +
    ked2_v2sigmalapl[0]*mgga_v3rhosigmatau[3] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[5] +
    mgga_v4rhosigmalapltau[7]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[9] + mgga_v4rhosigma2lapl[9];
  v4rhosigma2lapl[10] = ked1_v2rholapl[0]*mgga_v3sigma2tau[10] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigma2tau2[15] +
    mgga_v4sigma2lapltau[20]) + ked2_v2sigma2[0]*(ked1_v2rholapl[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[1]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked1_v2rholapl[0]*mgga_v3tau3[2] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[2]) + ked1_vlapl[0]*mgga_v4rhotau3[2] +
    mgga_v4rholapltau2[2]) + 2*ked2_vsigma[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[7] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[9] + mgga_v4sigmalapltau2[13]) + ked1_vlapl[0]*mgga_v4rhosigmatau2[7]
    + mgga_v4rhosigmalapltau[9]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[10] + mgga_v4rhosigma2lapl[10];
  v4rhosigma2lapl[11] = ked2_vlapl[0]*ked2_v2sigma2[0]*ked1_vrho[0]*mgga_v3tau3[2] +
    ked2_v2sigma2[0]*ked1_vrho[0]*mgga_v3lapltau2[4] + 2*ked2_v2sigmalapl[0]*ked1_vrho[0]*mgga_v3sigmatau2[7] +
    ked2_vlapl[0]*ked1_vrho[0]*mgga_v4sigma2tau2[16] + ked1_vrho[0]*mgga_v4sigma2lapltau[22] +
    ked2_v3sigma2lapl[0]*(ked1_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[1]) +
    ked2_vlapl[0]*ked2_v2sigma2[0]*mgga_v3rhotau2[2] + ked2_v2sigma2[0]*mgga_v3rholapltau[3] +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked1_vrho[0]*mgga_v4lapltau3[6] + ked2_vlapl[0]*(ked1_vrho[0]*mgga_v4tau4[3] +
    mgga_v4rhotau3[3]) + mgga_v4rholapltau2[5]) + 2*ked2_v2sigmalapl[0]*mgga_v3rhosigmatau[5] +
    2*ked2_vsigma[0]*(ked1_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[10] + mgga_v4sigmalapltau2[16]) +
    ked2_v2sigmalapl[0]*(ked1_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[2]) + ked2_vlapl[0]*mgga_v4rhosigmatau2[8] +
    mgga_v4rhosigmalapltau[11]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[11] + mgga_v4rhosigma2lapl[11];
  v4rhosigma2lapl[12] = ked1_vlapl[0]*ked1_v2sigma2[0]*ked2_vrho[0]*mgga_v3tau3[1] +
    ked1_v2sigma2[0]*ked2_vrho[0]*mgga_v3lapltau2[1] + 2*ked1_v2sigmalapl[0]*ked2_vrho[0]*mgga_v3sigmatau2[1] +
    ked1_vlapl[0]*ked2_vrho[0]*mgga_v4sigma2tau2[1] + ked2_vrho[0]*mgga_v4sigma2lapltau[1] +
    ked1_v3sigma2lapl[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    ked1_vlapl[0]*ked1_v2sigma2[0]*mgga_v3rhotau2[3] + ked1_v2sigma2[0]*mgga_v3rholapltau[4] +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked2_vrho[0]*mgga_v4lapltau3[1] + ked1_vlapl[0]*(ked2_vrho[0]*mgga_v4tau4[1] +
    mgga_v4rhotau3[4]) + mgga_v4rholapltau2[6]) + 2*ked1_v2sigmalapl[0]*mgga_v3rhosigmatau[6] +
    2*ked1_vsigma[0]*(ked2_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[1] + mgga_v4sigmalapltau2[1]) +
    ked1_v2sigmalapl[0]*(ked2_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[3]) + ked1_vlapl[0]*mgga_v4rhosigmatau2[9] +
    mgga_v4rhosigmalapltau[12]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[12] + mgga_v4rhosigma2lapl[12];
  v4rhosigma2lapl[13] = ked2_v2rholapl[0]*mgga_v3sigma2tau[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigma2tau2[2] +
    mgga_v4sigma2lapltau[3]) + ked1_v2sigma2[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3rhotau2[4] +
    mgga_v3rholapltau[6]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked2_v2rholapl[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4rhotau3[5] +
    mgga_v4rholapltau2[9]) + 2*ked1_vsigma[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[2] + mgga_v4sigmalapltau2[4]) + ked2_vlapl[0]*mgga_v4rhosigmatau2[10]
    + mgga_v4rhosigmalapltau[14]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[13] + mgga_v4rhosigma2lapl[13];
  v4rhosigma2lapl[14] = ked2_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[4] + ked1_vsigma[0]*mgga_v4sigmalapltau2[7]
    + ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[5] + mgga_v4sigma2tau2[4]) + mgga_v4sigma2lapltau[5]) +
    ked1_v2sigmalapl[0]*mgga_v3rhosigmatau[8] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[12] +
    mgga_v4rhosigmalapltau[16]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[14] + mgga_v4rhosigma2lapl[14];
  v4rhosigma2lapl[15] = ked2_v2rholapl[0]*mgga_v3sigma2tau[3] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigma2tau2[5] +
    mgga_v4sigma2lapltau[7]) + ked1_vsigma[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[6] + mgga_v4sigmalapltau2[10]) + ked2_vlapl[0]*mgga_v4rhosigmatau2[13]
    + mgga_v4rhosigmalapltau[18]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[15] + mgga_v4rhosigma2lapl[15];
  v4rhosigma2lapl[16] = ked2_v2rhosigma[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[1] + ked1_vsigma[0]*mgga_v3lapltau2[1] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + mgga_v3sigmalapltau[1]) +
    ked2_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[7] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[9] +
    mgga_v4sigmalapltau2[13]) + ked2_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[2] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[2]) + ked1_vlapl[0]*mgga_v4sigmatau3[2] +
    mgga_v4sigmalapltau2[2]) + ked1_vlapl[0]*mgga_v4sigma2tau2[7] + mgga_v4sigma2lapltau[9]) +
    ked1_v2sigmalapl[0]*(ked2_vsigma[0]*mgga_v3rhotau2[4] + mgga_v3rhosigmatau[10]) +
    ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v4rholapltau2[7] + ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4rhotau3[5] +
    mgga_v4rhosigmatau2[15]) + mgga_v4rhosigmalapltau[20]) + ked2_vsigma[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[10] +
    mgga_v4rhosigmalapltau[13]) + ked1_vlapl[0]*mgga_v4rhosigma2tau[16] + mgga_v4rhosigma2lapl[16];
  v4rhosigma2lapl[17] = ked2_v3rhosigmalapl[0]*(ked1_vsigma[0]*mgga_v2tau2[1] + mgga_v2sigmatau[1]) +
    ked2_v2rhosigma[0]*(ked1_vsigma[0]*mgga_v3lapltau2[4] + ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[2]) + mgga_v3sigmalapltau[3]) + ked2_v2rholapl[0]*(ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[7]) + ked2_vsigma[0]*mgga_v3sigmatau2[2] + mgga_v3sigma2tau[5]) +
    ked2_vrho[0]*(ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[3] +
    mgga_v4lapltau3[6]) + ked2_vlapl[0]*mgga_v4sigmatau3[10] + mgga_v4sigmalapltau2[16]) +
    ked2_v2sigmalapl[0]*mgga_v3sigmatau2[2] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[3] +
    mgga_v4sigmalapltau2[5]) + ked2_vlapl[0]*mgga_v4sigma2tau2[8] + mgga_v4sigma2lapltau[11]) +
    ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3rhotau2[4] + ked2_vsigma[0]*mgga_v4rholapltau2[10] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4rhotau3[6] + mgga_v4rhosigmatau2[16]) + mgga_v4rhosigmalapltau[22]) +
    ked2_v2sigmalapl[0]*mgga_v3rhosigmatau[7] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[11] +
    mgga_v4rhosigmalapltau[15]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[17] + mgga_v4rhosigma2lapl[17];
  v4rhosigma2lapl[18] = ked2_vrho[0]*mgga_v4sigma2lapltau[13] + ked1_vlapl[0]*(ked2_vrho[0]*mgga_v4sigma2tau2[10] +
    mgga_v4rhosigma2tau[18]) + mgga_v4rhosigma2lapl[18];
  v4rhosigma2lapl[19] = ked2_v2rholapl[0]*mgga_v3sigma2tau[7] + ked2_vrho[0]*mgga_v4sigma2lapltau[15] +
    ked2_vlapl[0]*(ked2_vrho[0]*mgga_v4sigma2tau2[11] + mgga_v4rhosigma2tau[19]) + mgga_v4rhosigma2lapl[19];
  v4rhosigma2lapl[20] = ked2_v2rhosigma[0]*mgga_v3sigmalapltau[5] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmalapltau2[8] + mgga_v4sigma2lapltau[17]) +
    ked2_vsigma[0]*mgga_v4rhosigmalapltau[17] + ked1_vlapl[0]*(ked2_v2rhosigma[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[6] + mgga_v4sigma2tau2[13]) + ked2_vsigma[0]*mgga_v4rhosigmatau2[13]
    + mgga_v4rhosigma2tau[20]) + mgga_v4rhosigma2lapl[20];
  v4rhosigma2lapl[21] = ked2_v3rhosigmalapl[0]*mgga_v2sigmatau[3] +
    ked2_v2rhosigma[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[5] + mgga_v3sigmalapltau[7]) +
    ked2_v2rholapl[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[5] + mgga_v3sigma2tau[9]) +
    ked2_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[5] + ked2_vsigma[0]*mgga_v4sigmalapltau2[11] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[7] + mgga_v4sigma2tau2[14]) + mgga_v4sigma2lapltau[19]) +
    ked2_v2sigmalapl[0]*mgga_v3rhosigmatau[9] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[14] +
    mgga_v4rhosigmalapltau[19]) + ked2_vlapl[0]*mgga_v4rhosigma2tau[21] + mgga_v4rhosigma2lapl[21];
  v4rhosigma2lapl[22] = ked2_v3rhosigma2[0]*mgga_v2lapltau[1] + 2*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3lapltau2[2]
    + ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4lapltau3[3] + 2*ked2_v2rhosigma[0]*mgga_v3sigmalapltau[9] +
    2*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigmalapltau2[14] + ked2_vrho[0]*mgga_v4sigma2lapltau[21] +
    ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3lapltau2[2] + mgga_v3rholapltau[5]) +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rholapltau2[8] + 2*ked2_vsigma[0]*mgga_v4rhosigmalapltau[21] +
    ked1_vlapl[0]*(ked2_v3rhosigma2[0]*mgga_v2tau2[1] + 2*ked2_v2rhosigma[0]*mgga_v3sigmatau2[7] +
    ked2_vrho[0]*mgga_v4sigma2tau2[16] + ked2_v2sigma2[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[6]) +
    2*ked2_vsigma[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[2] + ked2_vrho[0]*mgga_v4sigmatau3[10] + mgga_v4rhosigmatau2[16])
    + mgga_v4rhosigma2tau[22]) + mgga_v4rhosigma2lapl[22];
  v4rhosigma2lapl[23] = ked2_v4rhosigma2lapl[0]*mgga_vtau[1] + 2*ked2_v2sigmalapl[0]*ked2_v2rhosigma[0]*mgga_v2tau2[2]
    + 2*ked2_vsigma[0]*ked2_v3rhosigmalapl[0]*mgga_v2tau2[2] + ked2_vlapl[0]*ked2_v3rhosigma2[0]*mgga_v2tau2[2] +
    2*ked2_vsigma[0]*ked2_v2sigmalapl[0]*ked2_vrho[0]*mgga_v3tau3[3] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_v2rholapl[0]*mgga_v3tau3[3] +
    2*ked2_vlapl[0]*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3tau3[3] +
    ked2_vlapl[0]*ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4tau4[4] + ked2_v3rhosigma2[0]*mgga_v2lapltau[3] +
    2*ked2_vsigma[0]*ked2_v2rhosigma[0]*mgga_v3lapltau2[5] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4lapltau3[7] + 2*ked2_v3rhosigmalapl[0]*mgga_v2sigmatau[5] +
    2*ked2_v2sigmalapl[0]*ked2_vrho[0]*mgga_v3sigmatau2[8] + 2*ked2_vsigma[0]*ked2_v2rholapl[0]*mgga_v3sigmatau2[8] +
    2*ked2_vlapl[0]*ked2_v2rhosigma[0]*mgga_v3sigmatau2[8] +
    2*ked2_vlapl[0]*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigmatau3[11] + 2*ked2_v2rhosigma[0]*mgga_v3sigmalapltau[11] +
    2*ked2_vsigma[0]*ked2_vrho[0]*mgga_v4sigmalapltau2[17] + ked2_v2rholapl[0]*mgga_v3sigma2tau[11] +
    ked2_vlapl[0]*ked2_vrho[0]*mgga_v4sigma2tau2[17] + ked2_vrho[0]*mgga_v4sigma2lapltau[23] +
    ked2_v3sigma2lapl[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) +
    2*ked2_vsigma[0]*ked2_v2sigmalapl[0]*mgga_v3rhotau2[5] +
    ked2_vlapl[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rhotau3[7] + ked2_v2sigma2[0]*(ked2_v2rholapl[0]*mgga_v2tau2[2]
    + ked2_vrho[0]*mgga_v3lapltau2[5] + ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[3] + mgga_v3rhotau2[5]) +
    mgga_v3rholapltau[7]) + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4rholapltau2[11] +
    2*ked2_v2sigmalapl[0]*mgga_v3rhosigmatau[11] + 2*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4rhosigmatau2[17] +
    2*ked2_vsigma[0]*mgga_v4rhosigmalapltau[23] + ked2_vlapl[0]*mgga_v4rhosigma2tau[23] + mgga_v4rhosigma2lapl[23];
  v4rhosigma2tau[1] = 0;
  v4rhosigma2tau[2] = 0;
  v4rhosigma2tau[3] = 0;
  v4rhosigma2tau[4] = 0;
  v4rhosigma2tau[5] = 0;
  v4rhosigma2tau[6] = 0;
  v4rhosigma2tau[7] = 0;
  v4rhosigma2tau[8] = 0;
  v4rhosigma2tau[9] = 0;
  v4rhosigma2tau[10] = 0;
  v4rhosigma2tau[11] = 0;
  v4rhosigma2tau[12] = 0;
  v4rhosigma2tau[13] = 0;
  v4rhosigma2tau[14] = 0;
  v4rhosigma2tau[15] = 0;
  v4rhosigma2tau[16] = 0;
  v4rhosigma2tau[17] = 0;
  v4rhosigma2tau[18] = 0;
  v4rhosigma2tau[19] = 0;
  v4rhosigma2tau[20] = 0;
  v4rhosigma2tau[21] = 0;
  v4rhosigma2tau[22] = 0;
  v4rhosigma2tau[23] = 0;
  v4rhosigmalapl2[1] = ked1_v3rhosigmalapl[0]*(ked2_vlapl[0]*mgga_v2tau2[1] + mgga_v2lapltau[2]) +
    ked1_v2rhosigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[3]) +
    ked2_vlapl[0]*mgga_v3lapltau2[1] + mgga_v3lapl2tau[2]) + ked1_v2rholapl[0]*(ked1_vsigma[0]*mgga_v3lapltau2[3] +
    ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + mgga_v3sigmalapltau[2]) +
    ked1_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3lapltau2[3] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4lapltau3[4] +
    mgga_v4lapl2tau2[3]) + ked1_vlapl[0]*mgga_v4sigmalapltau2[3] + ked2_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) + ked1_vlapl[0]*mgga_v4sigmatau3[1] +
    mgga_v4sigmalapltau2[1]) + mgga_v4sigmalapl2tau[2]) + ked1_v2sigmalapl[0]*(ked2_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[2]) + ked1_vsigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhotau3[1] + mgga_v4rholapltau2[3]) +
    ked2_vlapl[0]*mgga_v4rholapltau2[1] + mgga_v4rholapl2tau[2]) + ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[1]
    + mgga_v4rhosigmalapltau[2]) + ked2_vlapl[0]*mgga_v4rhosigmalapltau[1] + mgga_v4rhosigmalapl2[1];
  v4rhosigmalapl2[2] = ked1_v2rhosigma[0]*mgga_v3lapl2tau[4] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4lapl2tau2[6] +
    mgga_v4sigmalapl2tau[4]) + ked1_vsigma[0]*mgga_v4rholapl2tau[4] +
    ked2_v2lapl2[0]*(ked1_v2rhosigma[0]*mgga_v2tau2[1] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v3tau3[1] +
    mgga_v3sigmatau2[1]) + ked1_vsigma[0]*mgga_v3rhotau2[1] + mgga_v3rhosigmatau[1]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked1_v2rhosigma[0]*mgga_v3tau3[2] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4tau4[2] +
    mgga_v4sigmatau3[2]) + ked1_vsigma[0]*mgga_v4rhotau3[2] + mgga_v4rhosigmatau2[2]) +
    2*ked2_vlapl[0]*(ked1_v2rhosigma[0]*mgga_v3lapltau2[4] + ked1_vrho[0]*(ked1_vsigma[0]*mgga_v4lapltau3[5] +
    mgga_v4sigmalapltau2[4]) + ked1_vsigma[0]*mgga_v4rholapltau2[4] + mgga_v4rhosigmalapltau[3]) +
    mgga_v4rhosigmalapl2[2];
  v4rhosigmalapl2[3] = ked1_v3rholapl2[0]*mgga_v2sigmatau[2] + 2*ked1_v2rholapl[0]*mgga_v3sigmalapltau[4] +
    ked1_vrho[0]*mgga_v4sigmalapl2tau[6] + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3sigmatau2[3] + mgga_v3rhosigmatau[2]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vrho[0]*mgga_v4sigmatau3[4] + mgga_v4rhosigmatau2[3]) +
    2*ked1_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[3] + ked1_vrho[0]*mgga_v4sigmalapltau2[6] +
    mgga_v4rhosigmalapltau[4]) + mgga_v4rhosigmalapl2[3];
  v4rhosigmalapl2[4] = ked1_v2rholapl[0]*mgga_v3sigmalapltau[6] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmalapltau2[9] +
    mgga_v4sigmalapl2tau[8]) + ked1_vlapl[0]*mgga_v4rhosigmalapltau[6] +
    ked2_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3sigmatau2[4] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[5] +
    mgga_v4sigmalapltau2[7]) + ked1_vlapl[0]*mgga_v4rhosigmatau2[4] + mgga_v4rhosigmalapltau[5]) +
    mgga_v4rhosigmalapl2[4];
  v4rhosigmalapl2[5] = ked1_vrho[0]*mgga_v4sigmalapl2tau[10] + ked2_v2lapl2[0]*(ked1_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[3]) + ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vrho[0]*mgga_v4sigmatau3[6] + mgga_v4rhosigmatau2[5]) +
    2*ked2_vlapl[0]*(ked1_vrho[0]*mgga_v4sigmalapltau2[10] + mgga_v4rhosigmalapltau[7]) + mgga_v4rhosigmalapl2[5];
  v4rhosigmalapl2[6] = ked1_v3rholapl2[0]*mgga_v2sigmatau[4] + 2*ked1_vlapl[0]*ked1_v2rholapl[0]*mgga_v3sigmatau2[6] +
    ked1_vlapl[0]*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4sigmatau3[8] + 2*ked1_v2rholapl[0]*mgga_v3sigmalapltau[8] +
    2*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4sigmalapltau2[12] + ked1_vrho[0]*mgga_v4sigmalapl2tau[12] +
    ked2_vsigma[0]*(ked1_v3rholapl2[0]*mgga_v2tau2[1] + 2*ked1_v2rholapl[0]*mgga_v3lapltau2[1] +
    ked1_vrho[0]*mgga_v4lapl2tau2[1] + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[1]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[1]) +
    2*ked1_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3tau3[1] + ked1_vrho[0]*mgga_v4lapltau3[1] + mgga_v4rholapltau2[1]) +
    mgga_v4rholapl2tau[1]) + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3sigmatau2[6] + mgga_v3rhosigmatau[4]) +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rhosigmatau2[6] + 2*ked1_vlapl[0]*mgga_v4rhosigmalapltau[8] +
    mgga_v4rhosigmalapl2[6];
  v4rhosigmalapl2[7] = ked1_v2rholapl[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] + ked2_vsigma[0]*mgga_v3lapltau2[4] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + mgga_v3sigmalapltau[10]) +
    ked1_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3lapltau2[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4lapltau3[2] +
    mgga_v4lapl2tau2[4]) + ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4sigmatau3[9] +
    mgga_v4sigmalapltau2[15]) + ked2_vlapl[0]*mgga_v4sigmalapltau2[13] + mgga_v4sigmalapl2tau[14]) +
    ked2_v2sigmalapl[0]*(ked1_vlapl[0]*mgga_v3rhotau2[1] + mgga_v3rholapltau[1]) +
    ked2_vsigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhotau3[2] + mgga_v4rholapltau2[4]) +
    ked2_vlapl[0]*mgga_v4rholapltau2[2] + mgga_v4rholapl2tau[3]) + ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[7]
    + mgga_v4rhosigmalapltau[10]) + ked2_vlapl[0]*mgga_v4rhosigmalapltau[9] + mgga_v4rhosigmalapl2[7];
  v4rhosigmalapl2[8] = ked1_vrho[0]*(ked2_v3sigmalapl2[0]*mgga_v2tau2[1] + 2*ked2_v2sigmalapl[0]*mgga_v3lapltau2[4] +
    ked2_vsigma[0]*mgga_v4lapl2tau2[7] + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) +
    2*ked2_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] + ked2_vsigma[0]*mgga_v4lapltau3[6] + mgga_v4sigmalapltau2[16])
    + mgga_v4sigmalapl2tau[16]) + ked2_v3sigmalapl2[0]*mgga_v2rhotau[1] +
    2*ked2_v2sigmalapl[0]*(ked2_vlapl[0]*mgga_v3rhotau2[2] + mgga_v3rholapltau[3]) +
    ked2_vsigma[0]*(ked2_v2lapl2[0]*mgga_v3rhotau2[2] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rhotau3[3] +
    2*ked2_vlapl[0]*mgga_v4rholapltau2[5] + mgga_v4rholapl2tau[5]) + ked2_v2lapl2[0]*mgga_v3rhosigmatau[5] +
    ked2_vlapl[0]*mgga_v4rhosigmalapltau[11] + ked2_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[8] +
    mgga_v4rhosigmalapltau[11]) + mgga_v4rhosigmalapl2[8];
  v4rhosigmalapl2[9] = ked2_vrho[0]*(ked1_v3sigmalapl2[0]*mgga_v2tau2[1] + 2*ked1_v2sigmalapl[0]*mgga_v3lapltau2[1] +
    ked1_vsigma[0]*mgga_v4lapl2tau2[1] + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) +
    2*ked1_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] + ked1_vsigma[0]*mgga_v4lapltau3[1] + mgga_v4sigmalapltau2[1])
    + mgga_v4sigmalapl2tau[1]) + ked1_v3sigmalapl2[0]*mgga_v2rhotau[2] +
    2*ked1_v2sigmalapl[0]*(ked1_vlapl[0]*mgga_v3rhotau2[3] + mgga_v3rholapltau[4]) +
    ked1_vsigma[0]*(ked1_v2lapl2[0]*mgga_v3rhotau2[3] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rhotau3[4] +
    2*ked1_vlapl[0]*mgga_v4rholapltau2[6] + mgga_v4rholapl2tau[6]) + ked1_v2lapl2[0]*mgga_v3rhosigmatau[6] +
    ked1_vlapl[0]*mgga_v4rhosigmalapltau[12] + ked1_vlapl[0]*(ked1_vlapl[0]*mgga_v4rhosigmatau2[9] +
    mgga_v4rhosigmalapltau[12]) + mgga_v4rhosigmalapl2[9];
  v4rhosigmalapl2[10] = ked2_v2rholapl[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[1] + ked1_vsigma[0]*mgga_v3lapltau2[1] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) + mgga_v3sigmalapltau[1]) +
    ked2_vrho[0]*(ked1_v2sigmalapl[0]*mgga_v3lapltau2[4] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4lapltau3[5] +
    mgga_v4lapl2tau2[4]) + ked1_vlapl[0]*mgga_v4sigmalapltau2[4] + ked2_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[2] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[2]) + ked1_vlapl[0]*mgga_v4sigmatau3[2] +
    mgga_v4sigmalapltau2[2]) + mgga_v4sigmalapl2tau[3]) + ked1_v2sigmalapl[0]*(ked2_vlapl[0]*mgga_v3rhotau2[4] +
    mgga_v3rholapltau[6]) + ked1_vsigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhotau3[5] + mgga_v4rholapltau2[9]) +
    ked2_vlapl[0]*mgga_v4rholapltau2[7] + mgga_v4rholapl2tau[8]) + ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[10]
    + mgga_v4rhosigmalapltau[14]) + ked2_vlapl[0]*mgga_v4rhosigmalapltau[13] + mgga_v4rhosigmalapl2[10];
  v4rhosigmalapl2[11] = ked2_v3rholapl2[0]*mgga_v2sigmatau[1] + 2*ked2_vlapl[0]*ked2_v2rholapl[0]*mgga_v3sigmatau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4sigmatau3[3] + 2*ked2_v2rholapl[0]*mgga_v3sigmalapltau[3] +
    2*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4sigmalapltau2[5] + ked2_vrho[0]*mgga_v4sigmalapl2tau[5] +
    ked1_vsigma[0]*(ked2_v3rholapl2[0]*mgga_v2tau2[1] + 2*ked2_v2rholapl[0]*mgga_v3lapltau2[4] +
    ked2_vrho[0]*mgga_v4lapl2tau2[7] + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[6]) +
    2*ked2_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3tau3[2] + ked2_vrho[0]*mgga_v4lapltau3[6] + mgga_v4rholapltau2[10]) +
    mgga_v4rholapl2tau[10]) + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3sigmatau2[2] + mgga_v3rhosigmatau[7]) +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rhosigmatau2[11] + 2*ked2_vlapl[0]*mgga_v4rhosigmalapltau[15] +
    mgga_v4rhosigmalapl2[11];
  v4rhosigmalapl2[12] = ked2_vrho[0]*mgga_v4sigmalapl2tau[7] + ked1_v2lapl2[0]*(ked2_vrho[0]*mgga_v3sigmatau2[4] +
    mgga_v3rhosigmatau[8]) + ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vrho[0]*mgga_v4sigmatau3[5] + mgga_v4rhosigmatau2[12]) +
    2*ked1_vlapl[0]*(ked2_vrho[0]*mgga_v4sigmalapltau2[7] + mgga_v4rhosigmalapltau[16]) + mgga_v4rhosigmalapl2[12];
  v4rhosigmalapl2[13] = ked2_v2rholapl[0]*mgga_v3sigmalapltau[5] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmalapltau2[8]
    + mgga_v4sigmalapl2tau[9]) + ked1_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[4] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[6] + mgga_v4sigmalapltau2[10]) + ked2_vlapl[0]*mgga_v4rhosigmatau2[13]
    + mgga_v4rhosigmalapltau[18]) + ked2_vlapl[0]*mgga_v4rhosigmalapltau[17] + mgga_v4rhosigmalapl2[13];
  v4rhosigmalapl2[14] = ked2_v3rholapl2[0]*mgga_v2sigmatau[3] + 2*ked2_v2rholapl[0]*mgga_v3sigmalapltau[7] +
    ked2_vrho[0]*mgga_v4sigmalapl2tau[11] + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3sigmatau2[5] + mgga_v3rhosigmatau[9])
    + ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vrho[0]*mgga_v4sigmatau3[7] + mgga_v4rhosigmatau2[14]) +
    2*ked2_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3sigmatau2[5] + ked2_vrho[0]*mgga_v4sigmalapltau2[11] +
    mgga_v4rhosigmalapltau[19]) + mgga_v4rhosigmalapl2[14];
  v4rhosigmalapl2[15] = ked2_v2rhosigma[0]*mgga_v3lapl2tau[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4lapl2tau2[2] +
    mgga_v4sigmalapl2tau[13]) + ked2_vsigma[0]*mgga_v4rholapl2tau[7] +
    ked1_v2lapl2[0]*(ked2_v2rhosigma[0]*mgga_v2tau2[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v3tau3[2] +
    mgga_v3sigmatau2[7]) + ked2_vsigma[0]*mgga_v3rhotau2[4] + mgga_v3rhosigmatau[10]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked2_v2rhosigma[0]*mgga_v3tau3[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4tau4[2] +
    mgga_v4sigmatau3[9]) + ked2_vsigma[0]*mgga_v4rhotau3[5] + mgga_v4rhosigmatau2[15]) +
    2*ked1_vlapl[0]*(ked2_v2rhosigma[0]*mgga_v3lapltau2[1] + ked2_vrho[0]*(ked2_vsigma[0]*mgga_v4lapltau3[2] +
    mgga_v4sigmalapltau2[13]) + ked2_vsigma[0]*mgga_v4rholapltau2[7] + mgga_v4rhosigmalapltau[20]) +
    mgga_v4rhosigmalapl2[15];
  v4rhosigmalapl2[16] = ked2_v3rhosigmalapl[0]*(ked1_vlapl[0]*mgga_v2tau2[1] + mgga_v2lapltau[1]) +
    ked2_v2rhosigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) +
    ked2_vlapl[0]*mgga_v3lapltau2[2] + mgga_v3lapl2tau[3]) + ked2_v2rholapl[0]*(ked2_vsigma[0]*mgga_v3lapltau2[2] +
    ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + mgga_v3sigmalapltau[9]) +
    ked2_vrho[0]*(ked2_v2sigmalapl[0]*mgga_v3lapltau2[2] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4lapltau3[3] +
    mgga_v4lapl2tau2[5]) + ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[3] + mgga_v4lapltau3[6]) + ked2_vlapl[0]*mgga_v4sigmatau3[10] +
    mgga_v4sigmalapltau2[16]) + ked2_vlapl[0]*mgga_v4sigmalapltau2[14] + mgga_v4sigmalapl2tau[15]) +
    ked2_v2sigmalapl[0]*(ked1_vlapl[0]*mgga_v3rhotau2[4] + mgga_v3rholapltau[5]) +
    ked2_vsigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhotau3[6] + mgga_v4rholapltau2[10]) +
    ked2_vlapl[0]*mgga_v4rholapltau2[8] + mgga_v4rholapl2tau[9]) + ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[16]
    + mgga_v4rhosigmalapltau[22]) + ked2_vlapl[0]*mgga_v4rhosigmalapltau[21] + mgga_v4rhosigmalapl2[16];
  v4rhosigmalapl2[17] = ked2_v4rhosigmalapl2[0]*mgga_vtau[1] + 2*ked2_v3rhosigmalapl[0]*(ked2_vlapl[0]*mgga_v2tau2[2] +
    mgga_v2lapltau[3]) + ked2_v2rhosigma[0]*(ked2_v2lapl2[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3tau3[3] + 2*ked2_vlapl[0]*mgga_v3lapltau2[5] + mgga_v3lapl2tau[5]) +
    ked2_v3rholapl2[0]*(ked2_vsigma[0]*mgga_v2tau2[2] + mgga_v2sigmatau[5]) +
    2*ked2_v2rholapl[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[2] + ked2_vsigma[0]*mgga_v3lapltau2[5] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[3] + mgga_v3sigmatau2[8]) + mgga_v3sigmalapltau[11]) +
    ked2_vrho[0]*(ked2_v3sigmalapl2[0]*mgga_v2tau2[2] + 2*ked2_v2sigmalapl[0]*mgga_v3lapltau2[5] +
    ked2_vsigma[0]*mgga_v4lapl2tau2[8] + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3tau3[3] + mgga_v3sigmatau2[8]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[4] + mgga_v4sigmatau3[11]) +
    2*ked2_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[3] + ked2_vsigma[0]*mgga_v4lapltau3[7] + mgga_v4sigmalapltau2[17])
    + mgga_v4sigmalapl2tau[17]) + ked2_v3sigmalapl2[0]*mgga_v2rhotau[3] +
    2*ked2_v2sigmalapl[0]*(ked2_vlapl[0]*mgga_v3rhotau2[5] + mgga_v3rholapltau[7]) +
    ked2_vsigma[0]*(ked2_v2lapl2[0]*mgga_v3rhotau2[5] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rhotau3[7] +
    2*ked2_vlapl[0]*mgga_v4rholapltau2[11] + mgga_v4rholapl2tau[11]) + ked2_v2lapl2[0]*mgga_v3rhosigmatau[11] +
    ked2_vlapl[0]*mgga_v4rhosigmalapltau[23] + ked2_vlapl[0]*(ked2_vlapl[0]*mgga_v4rhosigmatau2[17] +
    mgga_v4rhosigmalapltau[23]) + mgga_v4rhosigmalapl2[17];
  v4rhosigmalapltau[1] = 0;
  v4rhosigmalapltau[2] = 0;
  v4rhosigmalapltau[3] = 0;
  v4rhosigmalapltau[4] = 0;
  v4rhosigmalapltau[5] = 0;
  v4rhosigmalapltau[6] = 0;
  v4rhosigmalapltau[7] = 0;
  v4rhosigmalapltau[8] = 0;
  v4rhosigmalapltau[9] = 0;
  v4rhosigmalapltau[10] = 0;
  v4rhosigmalapltau[11] = 0;
  v4rhosigmalapltau[12] = 0;
  v4rhosigmalapltau[13] = 0;
  v4rhosigmalapltau[14] = 0;
  v4rhosigmalapltau[15] = 0;
  v4rhosigmalapltau[16] = 0;
  v4rhosigmalapltau[17] = 0;
  v4rhosigmalapltau[18] = 0;
  v4rhosigmalapltau[19] = 0;
  v4rhosigmalapltau[20] = 0;
  v4rhosigmalapltau[21] = 0;
  v4rhosigmalapltau[22] = 0;
  v4rhosigmalapltau[23] = 0;
  v4rhosigmatau2[1] = 0;
  v4rhosigmatau2[2] = 0;
  v4rhosigmatau2[3] = 0;
  v4rhosigmatau2[4] = 0;
  v4rhosigmatau2[5] = 0;
  v4rhosigmatau2[6] = 0;
  v4rhosigmatau2[7] = 0;
  v4rhosigmatau2[8] = 0;
  v4rhosigmatau2[9] = 0;
  v4rhosigmatau2[10] = 0;
  v4rhosigmatau2[11] = 0;
  v4rhosigmatau2[12] = 0;
  v4rhosigmatau2[13] = 0;
  v4rhosigmatau2[14] = 0;
  v4rhosigmatau2[15] = 0;
  v4rhosigmatau2[16] = 0;
  v4rhosigmatau2[17] = 0;
  v4rholapl3[1] = ked1_v3rholapl2[0]*mgga_v2lapltau[2] + 2*ked1_vlapl[0]*ked1_v2rholapl[0]*mgga_v3lapltau2[3] +
    ked1_vlapl[0]*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4lapltau3[4] + 2*ked1_v2rholapl[0]*mgga_v3lapl2tau[2] +
    2*ked1_vlapl[0]*ked1_vrho[0]*mgga_v4lapl2tau2[3] + ked1_vrho[0]*mgga_v4lapl3tau[2] +
    ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3lapltau2[3] + mgga_v3rholapltau[2]) +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4rholapltau2[3] + 2*ked1_vlapl[0]*mgga_v4rholapl2tau[2] +
    ked2_vlapl[0]*(ked1_v3rholapl2[0]*mgga_v2tau2[1] + 2*ked1_v2rholapl[0]*mgga_v3lapltau2[1] +
    ked1_vrho[0]*mgga_v4lapl2tau2[1] + ked1_v2lapl2[0]*(ked1_vrho[0]*mgga_v3tau3[1] + mgga_v3rhotau2[1]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[1]) +
    2*ked1_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3tau3[1] + ked1_vrho[0]*mgga_v4lapltau3[1] + mgga_v4rholapltau2[1]) +
    mgga_v4rholapl2tau[1]) + mgga_v4rholapl3[1];
  v4rholapl3[2] = ked1_v2rholapl[0]*mgga_v3lapl2tau[4] + ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4lapl2tau2[6] +
    mgga_v4lapl3tau[4]) + ked1_vlapl[0]*mgga_v4rholapl2tau[4] + ked2_v2lapl2[0]*(ked1_v2rholapl[0]*mgga_v2tau2[1] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3rhotau2[1] +
    mgga_v3rholapltau[1]) + ked2_vlapl[0]*ked2_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3tau3[2] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[2]) + ked1_vlapl[0]*mgga_v4rhotau3[2] +
    mgga_v4rholapltau2[2]) + 2*ked2_vlapl[0]*(ked1_v2rholapl[0]*mgga_v3lapltau2[4] +
    ked1_vrho[0]*(ked1_vlapl[0]*mgga_v4lapltau3[5] + mgga_v4lapl2tau2[4]) + ked1_vlapl[0]*mgga_v4rholapltau2[4] +
    mgga_v4rholapl2tau[3]) + mgga_v4rholapl3[2];
  v4rholapl3[3] = ked1_vrho[0]*(3*ked2_v2lapl2[0]*mgga_v3lapltau2[4] + mgga_v4lapl3tau[6]) +
    ked2_v3lapl3[0]*(ked1_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[1]) +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[3]) +
    3*ked2_v2lapl2[0]*mgga_v3rholapltau[3] + 3*ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vrho[0]*mgga_v4lapltau3[6] +
    mgga_v4rholapltau2[5]) + 3*ked2_vlapl[0]*(ked1_vrho[0]*(ked2_v2lapl2[0]*mgga_v3tau3[2] + mgga_v4lapl2tau2[7]) +
    ked2_v2lapl2[0]*mgga_v3rhotau2[2] + mgga_v4rholapl2tau[5]) + mgga_v4rholapl3[3];
  v4rholapl3[4] = ked2_vrho[0]*(3*ked1_v2lapl2[0]*mgga_v3lapltau2[1] + mgga_v4lapl3tau[1]) +
    ked1_v3lapl3[0]*(ked2_vrho[0]*mgga_v2tau2[1] + mgga_v2rhotau[2]) +
    ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vrho[0]*mgga_v4tau4[1] + mgga_v4rhotau3[4]) +
    3*ked1_v2lapl2[0]*mgga_v3rholapltau[4] + 3*ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vrho[0]*mgga_v4lapltau3[1] +
    mgga_v4rholapltau2[6]) + 3*ked1_vlapl[0]*(ked2_vrho[0]*(ked1_v2lapl2[0]*mgga_v3tau3[1] + mgga_v4lapl2tau2[1]) +
    ked1_v2lapl2[0]*mgga_v3rhotau2[3] + mgga_v4rholapl2tau[6]) + mgga_v4rholapl3[4];
  v4rholapl3[5] = ked2_v2rholapl[0]*mgga_v3lapl2tau[1] + ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4lapl2tau2[2] +
    mgga_v4lapl3tau[3]) + ked1_v2lapl2[0]*(ked2_v2rholapl[0]*mgga_v2tau2[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3rhotau2[4] +
    mgga_v3rholapltau[6]) + ked1_vlapl[0]*ked1_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3tau3[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4rhotau3[5] +
    mgga_v4rholapltau2[9]) + 2*ked1_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3lapltau2[1] +
    ked2_vrho[0]*(ked2_vlapl[0]*mgga_v4lapltau3[2] + mgga_v4lapl2tau2[4]) + ked2_vlapl[0]*mgga_v4rholapltau2[7] +
    mgga_v4rholapl2tau[8]) + ked2_vlapl[0]*mgga_v4rholapl2tau[7] + mgga_v4rholapl3[5];
  v4rholapl3[6] = ked2_v3rholapl2[0]*mgga_v2lapltau[1] + 2*ked2_vlapl[0]*ked2_v2rholapl[0]*mgga_v3lapltau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4lapltau3[3] + 2*ked2_v2rholapl[0]*mgga_v3lapl2tau[3] +
    2*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4lapl2tau2[5] + ked2_vrho[0]*mgga_v4lapl3tau[5] +
    ked1_vlapl[0]*(ked2_v3rholapl2[0]*mgga_v2tau2[1] + 2*ked2_v2rholapl[0]*mgga_v3lapltau2[4] +
    ked2_vrho[0]*mgga_v4lapl2tau2[7] + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3tau3[2] + mgga_v3rhotau2[4]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vrho[0]*mgga_v4tau4[3] + mgga_v4rhotau3[6]) +
    2*ked2_vlapl[0]*(ked2_v2rholapl[0]*mgga_v3tau3[2] + ked2_vrho[0]*mgga_v4lapltau3[6] + mgga_v4rholapltau2[10]) +
    mgga_v4rholapl2tau[10]) + ked2_v2lapl2[0]*(ked2_vrho[0]*mgga_v3lapltau2[2] + mgga_v3rholapltau[5]) +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rholapltau2[8] + 2*ked2_vlapl[0]*mgga_v4rholapl2tau[9] + mgga_v4rholapl3[6];
  v4rholapl3[7] = ked2_v4rholapl3[0]*mgga_vtau[1] + 3*ked2_vlapl[0]*ked2_v3rholapl2[0]*mgga_v2tau2[2] +
    3*ked2_vlapl[0]*ked2_vlapl[0]*ked2_v2rholapl[0]*mgga_v3tau3[3] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4tau4[4] + 3*ked2_v3rholapl2[0]*mgga_v2lapltau[3] +
    6*ked2_vlapl[0]*ked2_v2rholapl[0]*mgga_v3lapltau2[5] +
    3*ked2_vlapl[0]*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4lapltau3[7] + 3*ked2_v2rholapl[0]*mgga_v3lapl2tau[5] +
    3*ked2_vlapl[0]*ked2_vrho[0]*mgga_v4lapl2tau2[8] + ked2_vrho[0]*mgga_v4lapl3tau[7] +
    ked2_v3lapl3[0]*(ked2_vrho[0]*mgga_v2tau2[2] + mgga_v2rhotau[3]) +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rhotau3[7] + 3*ked2_v2lapl2[0]*(ked2_v2rholapl[0]*mgga_v2tau2[2] +
    ked2_vrho[0]*mgga_v3lapltau2[5] + ked2_vlapl[0]*(ked2_vrho[0]*mgga_v3tau3[3] + mgga_v3rhotau2[5]) +
    mgga_v3rholapltau[7]) + 3*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4rholapltau2[11] +
    3*ked2_vlapl[0]*mgga_v4rholapl2tau[11] + mgga_v4rholapl3[7];
  v4rholapl2tau[1] = 0;
  v4rholapl2tau[2] = 0;
  v4rholapl2tau[3] = 0;
  v4rholapl2tau[4] = 0;
  v4rholapl2tau[5] = 0;
  v4rholapl2tau[6] = 0;
  v4rholapl2tau[7] = 0;
  v4rholapl2tau[8] = 0;
  v4rholapl2tau[9] = 0;
  v4rholapl2tau[10] = 0;
  v4rholapl2tau[11] = 0;
  v4rholapltau2[1] = 0;
  v4rholapltau2[2] = 0;
  v4rholapltau2[3] = 0;
  v4rholapltau2[4] = 0;
  v4rholapltau2[5] = 0;
  v4rholapltau2[6] = 0;
  v4rholapltau2[7] = 0;
  v4rholapltau2[8] = 0;
  v4rholapltau2[9] = 0;
  v4rholapltau2[10] = 0;
  v4rholapltau2[11] = 0;
  v4rhotau3[1] = 0;
  v4rhotau3[2] = 0;
  v4rhotau3[3] = 0;
  v4rhotau3[4] = 0;
  v4rhotau3[5] = 0;
  v4rhotau3[6] = 0;
  v4rhotau3[7] = 0;
  v4sigma4[1] = ked1_v3sigma3[0]*mgga_v2sigmatau[2] + ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[4]
    + 3*ked1_v2sigma2[0]*mgga_v3sigma2tau[2] + 3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigma2tau2[3] +
    3*ked1_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3sigmatau2[3] + mgga_v4sigma3tau[2]) + mgga_v4sigma4[1];
  v4sigma4[2] = ked1_v3sigma3[0]*mgga_v2sigmatau[4] + ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[8]
    + 3*ked1_v2sigma2[0]*mgga_v3sigma2tau[4] + 3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigma2tau2[6] +
    3*ked1_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3sigmatau2[6] + mgga_v4sigma3tau[4]) +
    ked2_vsigma[0]*(ked1_v3sigma3[0]*mgga_v2tau2[1] + ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4tau4[1] +
    3*ked1_v2sigma2[0]*mgga_v3sigmatau2[1] + 3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[1] +
    3*ked1_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3tau3[1] + mgga_v4sigma2tau2[1]) + mgga_v4sigma3tau[1]) +
    mgga_v4sigma4[2];
  v4sigma4[3] = ked1_v2sigma2[0]*mgga_v3sigma2tau[6] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigma2tau2[9] +
    2*ked1_vsigma[0]*mgga_v4sigma3tau[6] + mgga_v4sigma4[3];
  v4sigma4[4] = ked1_v2sigma2[0]*mgga_v3sigma2tau[8] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigma2tau2[12] +
    2*ked1_vsigma[0]*mgga_v4sigma3tau[8] + ked2_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3sigmatau2[4] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[5] + 2*ked1_vsigma[0]*mgga_v4sigma2tau2[4] + mgga_v4sigma3tau[3]) +
    mgga_v4sigma4[4];
  v4sigma4[5] = ked1_v2sigma2[0]*(ked2_v2sigma2[0]*mgga_v2tau2[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[2] +
    2*ked2_vsigma[0]*mgga_v3sigmatau2[7] + mgga_v3sigma2tau[10]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3tau3[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[2] +
    2*ked2_vsigma[0]*mgga_v4sigmatau3[9] + mgga_v4sigma2tau2[15]) +
    2*ked1_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3sigmatau2[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[2] +
    2*ked2_vsigma[0]*mgga_v4sigma2tau2[7] + mgga_v4sigma3tau[10]) + ked2_v2sigma2[0]*mgga_v3sigma2tau[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigma2tau2[2] + 2*ked2_vsigma[0]*mgga_v4sigma3tau[5] + mgga_v4sigma4[5];
  v4sigma4[6] = ked1_vsigma[0]*mgga_v4sigma3tau[12] + mgga_v4sigma4[6];
  v4sigma4[7] = ked1_vsigma[0]*(ked2_vsigma[0]*mgga_v4sigma2tau2[10] + mgga_v4sigma3tau[14]) +
    ked2_vsigma[0]*mgga_v4sigma3tau[7] + mgga_v4sigma4[7];
  v4sigma4[8] = ked1_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3sigmatau2[4] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[6] + 2*ked2_vsigma[0]*mgga_v4sigma2tau2[13] + mgga_v4sigma3tau[16])
    + ked2_v2sigma2[0]*mgga_v3sigma2tau[3] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigma2tau2[5] +
    2*ked2_vsigma[0]*mgga_v4sigma3tau[9] + mgga_v4sigma4[8];
  v4sigma4[9] = ked1_vsigma[0]*(ked2_v3sigma3[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[3] + 3*ked2_v2sigma2[0]*mgga_v3sigmatau2[7] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[10] + 3*ked2_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3tau3[2] +
    mgga_v4sigma2tau2[16]) + mgga_v4sigma3tau[18]) + ked2_v3sigma3[0]*mgga_v2sigmatau[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[3] + 3*ked2_v2sigma2[0]*mgga_v3sigma2tau[5] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigma2tau2[8] + 3*ked2_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3sigmatau2[2] +
    mgga_v4sigma3tau[11]) + mgga_v4sigma4[9];
  v4sigma4[10] = mgga_v4sigma4[10];
  v4sigma4[11] = ked2_vsigma[0]*mgga_v4sigma3tau[13] + mgga_v4sigma4[11];
  v4sigma4[12] = ked2_v2sigma2[0]*mgga_v3sigma2tau[7] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigma2tau2[11] +
    2*ked2_vsigma[0]*mgga_v4sigma3tau[15] + mgga_v4sigma4[12];
  v4sigma4[13] = ked2_v3sigma3[0]*mgga_v2sigmatau[3] + ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[7]
    + 3*ked2_v2sigma2[0]*mgga_v3sigma2tau[9] + 3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigma2tau2[14] +
    3*ked2_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3sigmatau2[5] + mgga_v4sigma3tau[17]) + mgga_v4sigma4[13];
  v4sigma4[14] = ked2_v4sigma4[0]*mgga_vtau[1] + 3*ked2_v2sigma2[0]*ked2_v2sigma2[0]*mgga_v2tau2[2] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[4] + 4*ked2_v3sigma3[0]*mgga_v2sigmatau[5]
    + 4*ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[11] +
    6*ked2_v2sigma2[0]*(ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[3] + 2*ked2_vsigma[0]*mgga_v3sigmatau2[8] +
    mgga_v3sigma2tau[11]) + 6*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigma2tau2[17] +
    4*ked2_vsigma[0]*(ked2_v3sigma3[0]*mgga_v2tau2[2] + mgga_v4sigma3tau[19]) + mgga_v4sigma4[14];
  v4sigma3lapl[1] = ked1_v3sigma3[0]*mgga_v2lapltau[2] +
    ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4lapltau3[4] + 3*ked1_v2sigma2[0]*mgga_v3sigmalapltau[2] +
    3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmalapltau2[3] + 3*ked1_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3lapltau2[3] +
    mgga_v4sigma2lapltau[2]) + ked2_vlapl[0]*(ked1_v3sigma3[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4tau4[1] + 3*ked1_v2sigma2[0]*mgga_v3sigmatau2[1] +
    3*ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[1] + 3*ked1_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3tau3[1] +
    mgga_v4sigma2tau2[1]) + mgga_v4sigma3tau[1]) + mgga_v4sigma3lapl[1];
  v4sigma3lapl[2] = ked1_v3sigma2lapl[0]*mgga_v2sigmatau[2] + ked1_v2sigma2[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[3] +
    mgga_v3sigmalapltau[4]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[4] +
    mgga_v4sigmalapltau2[6]) + 2*ked1_v2sigmalapl[0]*mgga_v3sigma2tau[2] +
    2*ked1_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[3] + ked1_vlapl[0]*mgga_v4sigma2tau2[3] +
    mgga_v4sigma2lapltau[4]) + ked1_vlapl[0]*mgga_v4sigma3tau[2] + mgga_v4sigma3lapl[2];
  v4sigma3lapl[3] = ked1_v2sigma2[0]*mgga_v3sigmalapltau[6] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmalapltau2[9] +
    2*ked1_vsigma[0]*mgga_v4sigma2lapltau[6] + ked2_vlapl[0]*(ked1_v2sigma2[0]*mgga_v3sigmatau2[4] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[5] + 2*ked1_vsigma[0]*mgga_v4sigma2tau2[4] + mgga_v4sigma3tau[3]) +
    mgga_v4sigma3lapl[3];
  v4sigma3lapl[4] = ked1_v3sigma2lapl[0]*mgga_v2sigmatau[4] + ked1_v2sigma2[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[6] +
    mgga_v3sigmalapltau[8]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[8] +
    mgga_v4sigmalapltau2[12]) + 2*ked1_v2sigmalapl[0]*mgga_v3sigma2tau[4] +
    2*ked1_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[6] + ked1_vlapl[0]*mgga_v4sigma2tau2[6] +
    mgga_v4sigma2lapltau[8]) + ked2_vsigma[0]*(ked1_v3sigma2lapl[0]*mgga_v2tau2[1] +
    ked1_v2sigma2[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) +
    ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) +
    2*ked1_v2sigmalapl[0]*mgga_v3sigmatau2[1] + 2*ked1_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked1_vlapl[0]*mgga_v4sigmatau3[1] + mgga_v4sigmalapltau2[1]) + ked1_vlapl[0]*mgga_v4sigma2tau2[1] +
    mgga_v4sigma2lapltau[1]) + ked1_vlapl[0]*mgga_v4sigma3tau[4] + mgga_v4sigma3lapl[4];
  v4sigma3lapl[5] = ked2_vsigma[0]*ked1_v2sigma2[0]*mgga_v3lapltau2[4] +
    ked1_vsigma[0]*ked1_vsigma[0]*ked2_vsigma[0]*mgga_v4lapltau3[5] + ked1_v2sigma2[0]*mgga_v3sigmalapltau[10] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmalapltau2[15] + 2*ked1_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmalapltau2[4] +
    2*ked1_vsigma[0]*mgga_v4sigma2lapltau[10] + ked2_v2sigmalapl[0]*(ked1_v2sigma2[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[1] + 2*ked1_vsigma[0]*mgga_v3sigmatau2[1] + mgga_v3sigma2tau[1]) +
    ked2_vsigma[0]*mgga_v4sigma2lapltau[3] + ked2_vlapl[0]*(ked1_v2sigma2[0]*mgga_v3sigmatau2[7] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4sigmatau3[9] + 2*ked1_vsigma[0]*mgga_v4sigma2tau2[7] +
    ked2_vsigma[0]*(ked1_v2sigma2[0]*mgga_v3tau3[2] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4tau4[2] +
    2*ked1_vsigma[0]*mgga_v4sigmatau3[2] + mgga_v4sigma2tau2[2]) + mgga_v4sigma3tau[5]) + mgga_v4sigma3lapl[5];
  v4sigma3lapl[6] = ked1_v2sigmalapl[0]*mgga_v3sigma2tau[6] + ked1_vsigma[0]*mgga_v4sigma2lapltau[12] +
    ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigma2tau2[9] + mgga_v4sigma3tau[6]) + mgga_v4sigma3lapl[6];
  v4sigma3lapl[7] = ked1_vsigma[0]*mgga_v4sigma2lapltau[14] + ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigma2tau2[10] +
    mgga_v4sigma3tau[7]) + mgga_v4sigma3lapl[7];
  v4sigma3lapl[8] = ked1_v2sigmalapl[0]*mgga_v3sigma2tau[8] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigma2tau2[12] +
    mgga_v4sigma2lapltau[16]) + ked2_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[4] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[5] + mgga_v4sigmalapltau2[7]) + ked1_vlapl[0]*mgga_v4sigma2tau2[4] +
    mgga_v4sigma2lapltau[5]) + ked1_vlapl[0]*mgga_v4sigma3tau[8] + mgga_v4sigma3lapl[8];
  v4sigma3lapl[9] = ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[4] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[6] + mgga_v4sigmalapltau2[10]) + ked2_vlapl[0]*mgga_v4sigma2tau2[13]
    + mgga_v4sigma2lapltau[18]) + ked2_v2sigmalapl[0]*mgga_v3sigma2tau[3] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigma2tau2[5] + mgga_v4sigma2lapltau[7]) + ked2_vlapl[0]*mgga_v4sigma3tau[9] +
    mgga_v4sigma3lapl[9];
  v4sigma3lapl[10] = ked1_vsigma[0]*ked2_v2sigma2[0]*mgga_v3lapltau2[1] +
    ked1_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4lapltau3[2] +
    2*ked1_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmalapltau2[13] + ked1_v2sigmalapl[0]*(ked2_v2sigma2[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[2] + 2*ked2_vsigma[0]*mgga_v3sigmatau2[7] + mgga_v3sigma2tau[10]) +
    ked1_vsigma[0]*mgga_v4sigma2lapltau[20] + ked2_v2sigma2[0]*mgga_v3sigmalapltau[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmalapltau2[2] + 2*ked2_vsigma[0]*mgga_v4sigma2lapltau[9] +
    ked1_vlapl[0]*(ked1_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3tau3[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[2] +
    2*ked2_vsigma[0]*mgga_v4sigmatau3[9] + mgga_v4sigma2tau2[15]) + ked2_v2sigma2[0]*mgga_v3sigmatau2[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[2] + 2*ked2_vsigma[0]*mgga_v4sigma2tau2[7] + mgga_v4sigma3tau[10]) +
    mgga_v4sigma3lapl[10];
  v4sigma3lapl[11] = ked1_vsigma[0]*(ked2_v3sigma2lapl[0]*mgga_v2tau2[1] +
    ked2_v2sigma2[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[3] + mgga_v4lapltau3[6]) +
    2*ked2_v2sigmalapl[0]*mgga_v3sigmatau2[7] + 2*ked2_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] +
    ked2_vlapl[0]*mgga_v4sigmatau3[10] + mgga_v4sigmalapltau2[16]) + ked2_vlapl[0]*mgga_v4sigma2tau2[16] +
    mgga_v4sigma2lapltau[22]) + ked2_v3sigma2lapl[0]*mgga_v2sigmatau[1] +
    ked2_v2sigma2[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[2] + mgga_v3sigmalapltau[3]) +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[3] + mgga_v4sigmalapltau2[5]) +
    2*ked2_v2sigmalapl[0]*mgga_v3sigma2tau[5] + 2*ked2_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[2] +
    ked2_vlapl[0]*mgga_v4sigma2tau2[8] + mgga_v4sigma2lapltau[11]) + ked2_vlapl[0]*mgga_v4sigma3tau[11] +
    mgga_v4sigma3lapl[11];
  v4sigma3lapl[12] = ked1_vlapl[0]*mgga_v4sigma3tau[12] + mgga_v4sigma3lapl[12];
  v4sigma3lapl[13] = ked2_vlapl[0]*mgga_v4sigma3tau[13] + mgga_v4sigma3lapl[13];
  v4sigma3lapl[14] = ked2_vsigma[0]*mgga_v4sigma2lapltau[13] + ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigma2tau2[10] +
    mgga_v4sigma3tau[14]) + mgga_v4sigma3lapl[14];
  v4sigma3lapl[15] = ked2_v2sigmalapl[0]*mgga_v3sigma2tau[7] + ked2_vsigma[0]*mgga_v4sigma2lapltau[15] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigma2tau2[11] + mgga_v4sigma3tau[15]) + mgga_v4sigma3lapl[15];
  v4sigma3lapl[16] = ked2_v2sigma2[0]*mgga_v3sigmalapltau[5] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmalapltau2[8] +
    2*ked2_vsigma[0]*mgga_v4sigma2lapltau[17] + ked1_vlapl[0]*(ked2_v2sigma2[0]*mgga_v3sigmatau2[4] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[6] + 2*ked2_vsigma[0]*mgga_v4sigma2tau2[13] + mgga_v4sigma3tau[16])
    + mgga_v4sigma3lapl[16];
  v4sigma3lapl[17] = ked2_v3sigma2lapl[0]*mgga_v2sigmatau[3] + ked2_v2sigma2[0]*(ked2_vlapl[0]*mgga_v3sigmatau2[5] +
    mgga_v3sigmalapltau[7]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[7] +
    mgga_v4sigmalapltau2[11]) + 2*ked2_v2sigmalapl[0]*mgga_v3sigma2tau[9] +
    2*ked2_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[5] + ked2_vlapl[0]*mgga_v4sigma2tau2[14] +
    mgga_v4sigma2lapltau[19]) + ked2_vlapl[0]*mgga_v4sigma3tau[17] + mgga_v4sigma3lapl[17];
  v4sigma3lapl[18] = ked2_v3sigma3[0]*mgga_v2lapltau[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4lapltau3[3] + 3*ked2_v2sigma2[0]*mgga_v3sigmalapltau[9] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmalapltau2[14] + 3*ked2_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3lapltau2[2] +
    mgga_v4sigma2lapltau[21]) + ked1_vlapl[0]*(ked2_v3sigma3[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[3] + 3*ked2_v2sigma2[0]*mgga_v3sigmatau2[7] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4sigmatau3[10] + 3*ked2_vsigma[0]*(ked2_v2sigma2[0]*mgga_v3tau3[2] +
    mgga_v4sigma2tau2[16]) + mgga_v4sigma3tau[18]) + mgga_v4sigma3lapl[18];
  v4sigma3lapl[19] = ked2_v4sigma3lapl[0]*mgga_vtau[1] + ked2_vlapl[0]*ked2_v3sigma3[0]*mgga_v2tau2[2] +
    ked2_v3sigma3[0]*mgga_v2lapltau[3] + ked2_vsigma[0]*ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[4] +
    mgga_v4lapltau3[7]) + 3*ked2_v3sigma2lapl[0]*mgga_v2sigmatau[5] +
    3*ked2_vlapl[0]*ked2_v2sigma2[0]*mgga_v3sigmatau2[8] + 3*ked2_v2sigma2[0]*mgga_v3sigmalapltau[11] +
    3*ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[11] + mgga_v4sigmalapltau2[17]) +
    3*ked2_v2sigmalapl[0]*(ked2_v2sigma2[0]*mgga_v2tau2[2] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[3] +
    2*ked2_vsigma[0]*mgga_v3sigmatau2[8] + mgga_v3sigma2tau[11]) +
    3*ked2_vsigma[0]*(ked2_v3sigma2lapl[0]*mgga_v2tau2[2] + ked2_v2sigma2[0]*mgga_v3lapltau2[5] +
    ked2_vlapl[0]*(ked2_v2sigma2[0]*mgga_v3tau3[3] + mgga_v4sigma2tau2[17]) + mgga_v4sigma2lapltau[23]) +
    ked2_vlapl[0]*mgga_v4sigma3tau[19] + mgga_v4sigma3lapl[19];
  v4sigma3tau[1] = 0;
  v4sigma3tau[2] = 0;
  v4sigma3tau[3] = 0;
  v4sigma3tau[4] = 0;
  v4sigma3tau[5] = 0;
  v4sigma3tau[6] = 0;
  v4sigma3tau[7] = 0;
  v4sigma3tau[8] = 0;
  v4sigma3tau[9] = 0;
  v4sigma3tau[10] = 0;
  v4sigma3tau[11] = 0;
  v4sigma3tau[12] = 0;
  v4sigma3tau[13] = 0;
  v4sigma3tau[14] = 0;
  v4sigma3tau[15] = 0;
  v4sigma3tau[16] = 0;
  v4sigma3tau[17] = 0;
  v4sigma3tau[18] = 0;
  v4sigma3tau[19] = 0;
  v4sigma2lapl2[1] = ked1_v3sigma2lapl[0]*mgga_v2lapltau[2] + ked1_v2sigma2[0]*(ked1_vlapl[0]*mgga_v3lapltau2[3] +
    mgga_v3lapl2tau[2]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4lapltau3[4] + mgga_v4lapl2tau2[3]) +
    2*ked1_v2sigmalapl[0]*mgga_v3sigmalapltau[2] + 2*ked1_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3lapltau2[3] +
    ked1_vlapl[0]*mgga_v4sigmalapltau2[3] + mgga_v4sigmalapl2tau[2]) + ked1_vlapl[0]*mgga_v4sigma2lapltau[2] +
    ked2_vlapl[0]*(ked1_v3sigma2lapl[0]*mgga_v2tau2[1] + ked1_v2sigma2[0]*(ked1_vlapl[0]*mgga_v3tau3[1] +
    mgga_v3lapltau2[1]) + ked1_vsigma[0]*ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[1] + mgga_v4lapltau3[1]) +
    2*ked1_v2sigmalapl[0]*mgga_v3sigmatau2[1] + 2*ked1_vsigma[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked1_vlapl[0]*mgga_v4sigmatau3[1] + mgga_v4sigmalapltau2[1]) + ked1_vlapl[0]*mgga_v4sigma2tau2[1] +
    mgga_v4sigma2lapltau[1]) + mgga_v4sigma2lapl2[1];
  v4sigma2lapl2[2] = ked1_v2sigma2[0]*mgga_v3lapl2tau[4] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4lapl2tau2[6] +
    2*ked1_vsigma[0]*mgga_v4sigmalapl2tau[4] + ked2_v2lapl2[0]*(ked1_v2sigma2[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*ked1_vsigma[0]*mgga_v3tau3[1] + 2*ked1_vsigma[0]*mgga_v3sigmatau2[1] + mgga_v3sigma2tau[1]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked1_v2sigma2[0]*mgga_v3tau3[2] + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4tau4[2] +
    2*ked1_vsigma[0]*mgga_v4sigmatau3[2] + mgga_v4sigma2tau2[2]) + 2*ked2_vlapl[0]*(ked1_v2sigma2[0]*mgga_v3lapltau2[4]
    + ked1_vsigma[0]*ked1_vsigma[0]*mgga_v4lapltau3[5] + 2*ked1_vsigma[0]*mgga_v4sigmalapltau2[4] +
    mgga_v4sigma2lapltau[3]) + mgga_v4sigma2lapl2[2];
  v4sigma2lapl2[3] = ked1_v3sigmalapl2[0]*mgga_v2sigmatau[2] + 2*ked1_v2sigmalapl[0]*mgga_v3sigmalapltau[4] +
    ked1_vsigma[0]*mgga_v4sigmalapl2tau[6] + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[3] + mgga_v3sigma2tau[2])
    + ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[4] + mgga_v4sigma2tau2[3]) +
    2*ked1_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[3] + ked1_vsigma[0]*mgga_v4sigmalapltau2[6] +
    mgga_v4sigma2lapltau[4]) + mgga_v4sigma2lapl2[3];
  v4sigma2lapl2[4] = ked1_v2sigmalapl[0]*mgga_v3sigmalapltau[6] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmalapltau2[9]
    + mgga_v4sigmalapl2tau[8]) + ked1_vlapl[0]*mgga_v4sigma2lapltau[6] +
    ked2_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3sigmatau2[4] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4sigmatau3[5] +
    mgga_v4sigmalapltau2[7]) + ked1_vlapl[0]*mgga_v4sigma2tau2[4] + mgga_v4sigma2lapltau[5]) + mgga_v4sigma2lapl2[4];
  v4sigma2lapl2[5] = ked1_vsigma[0]*mgga_v4sigmalapl2tau[10] + ked2_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[4] +
    mgga_v3sigma2tau[3]) + ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigmatau3[6] + mgga_v4sigma2tau2[5]) +
    2*ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v4sigmalapltau2[10] + mgga_v4sigma2lapltau[7]) + mgga_v4sigma2lapl2[5];
  v4sigma2lapl2[6] = ked1_v3sigmalapl2[0]*mgga_v2sigmatau[4] + 2*ked1_vlapl[0]*ked1_v2sigmalapl[0]*mgga_v3sigmatau2[6]
    + ked1_vlapl[0]*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4sigmatau3[8] + 2*ked1_v2sigmalapl[0]*mgga_v3sigmalapltau[8] +
    2*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4sigmalapltau2[12] + ked1_vsigma[0]*mgga_v4sigmalapl2tau[12] +
    ked2_vsigma[0]*(ked1_v3sigmalapl2[0]*mgga_v2tau2[1] + 2*ked1_v2sigmalapl[0]*mgga_v3lapltau2[1] +
    ked1_vsigma[0]*mgga_v4lapl2tau2[1] + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) +
    2*ked1_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] + ked1_vsigma[0]*mgga_v4lapltau3[1] + mgga_v4sigmalapltau2[1])
    + mgga_v4sigmalapl2tau[1]) + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3sigmatau2[6] + mgga_v3sigma2tau[4]) +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigma2tau2[6] + 2*ked1_vlapl[0]*mgga_v4sigma2lapltau[8] + mgga_v4sigma2lapl2[6];
  v4sigma2lapl2[7] = ked1_v2sigmalapl[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] + ked2_vsigma[0]*mgga_v3lapltau2[4] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) + mgga_v3sigmalapltau[10]) +
    ked1_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3lapltau2[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4lapltau3[2] +
    mgga_v4lapl2tau2[4]) + ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4sigmatau3[9] +
    mgga_v4sigmalapltau2[15]) + ked2_vlapl[0]*mgga_v4sigmalapltau2[13] + mgga_v4sigmalapl2tau[14]) +
    ked2_v2sigmalapl[0]*(ked1_vlapl[0]*mgga_v3sigmatau2[1] + mgga_v3sigmalapltau[1]) +
    ked2_vsigma[0]*(ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[2] + mgga_v4sigmalapltau2[4]) +
    ked2_vlapl[0]*mgga_v4sigmalapltau2[2] + mgga_v4sigmalapl2tau[3]) +
    ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4sigma2tau2[7] + mgga_v4sigma2lapltau[10]) +
    ked2_vlapl[0]*mgga_v4sigma2lapltau[9] + mgga_v4sigma2lapl2[7];
  v4sigma2lapl2[8] = ked1_vsigma[0]*(ked2_v3sigmalapl2[0]*mgga_v2tau2[1] + 2*ked2_v2sigmalapl[0]*mgga_v3lapltau2[4] +
    ked2_vsigma[0]*mgga_v4lapl2tau2[7] + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) +
    2*ked2_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] + ked2_vsigma[0]*mgga_v4lapltau3[6] + mgga_v4sigmalapltau2[16])
    + mgga_v4sigmalapl2tau[16]) + ked2_v3sigmalapl2[0]*mgga_v2sigmatau[1] +
    2*ked2_vlapl[0]*ked2_v2sigmalapl[0]*mgga_v3sigmatau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4sigmatau3[3] + 2*ked2_v2sigmalapl[0]*mgga_v3sigmalapltau[3] +
    2*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4sigmalapltau2[5] + ked2_vsigma[0]*mgga_v4sigmalapl2tau[5] +
    ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[2] + mgga_v3sigma2tau[5]) +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigma2tau2[8] + 2*ked2_vlapl[0]*mgga_v4sigma2lapltau[11] +
    mgga_v4sigma2lapl2[8];
  v4sigma2lapl2[9] = ked1_v2lapl2[0]*mgga_v3sigma2tau[6] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigma2tau2[9] +
    2*ked1_vlapl[0]*mgga_v4sigma2lapltau[12] + mgga_v4sigma2lapl2[9];
  v4sigma2lapl2[10] = ked1_vlapl[0]*(ked2_vlapl[0]*mgga_v4sigma2tau2[10] + mgga_v4sigma2lapltau[14]) +
    ked2_vlapl[0]*mgga_v4sigma2lapltau[13] + mgga_v4sigma2lapl2[10];
  v4sigma2lapl2[11] = ked2_v2lapl2[0]*mgga_v3sigma2tau[7] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigma2tau2[11] +
    2*ked2_vlapl[0]*mgga_v4sigma2lapltau[15] + mgga_v4sigma2lapl2[11];
  v4sigma2lapl2[12] = ked2_vsigma[0]*mgga_v4sigmalapl2tau[7] + ked1_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[4] +
    mgga_v3sigma2tau[8]) + ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[5] + mgga_v4sigma2tau2[12]) +
    2*ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigmalapltau2[7] + mgga_v4sigma2lapltau[16]) + mgga_v4sigma2lapl2[12];
  v4sigma2lapl2[13] = ked2_v2sigmalapl[0]*mgga_v3sigmalapltau[5] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmalapltau2[8] + mgga_v4sigmalapl2tau[9]) +
    ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[4] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4sigmatau3[6] +
    mgga_v4sigmalapltau2[10]) + ked2_vlapl[0]*mgga_v4sigma2tau2[13] + mgga_v4sigma2lapltau[18]) +
    ked2_vlapl[0]*mgga_v4sigma2lapltau[17] + mgga_v4sigma2lapl2[13];
  v4sigma2lapl2[14] = ked2_v3sigmalapl2[0]*mgga_v2sigmatau[3] + 2*ked2_v2sigmalapl[0]*mgga_v3sigmalapltau[7] +
    ked2_vsigma[0]*mgga_v4sigmalapl2tau[11] + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3sigmatau2[5] +
    mgga_v3sigma2tau[9]) + ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4sigmatau3[7] + mgga_v4sigma2tau2[14]) +
    2*ked2_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3sigmatau2[5] + ked2_vsigma[0]*mgga_v4sigmalapltau2[11] +
    mgga_v4sigma2lapltau[19]) + mgga_v4sigma2lapl2[14];
  v4sigma2lapl2[15] = ked2_v2sigma2[0]*mgga_v3lapl2tau[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4lapl2tau2[2] +
    2*ked2_vsigma[0]*mgga_v4sigmalapl2tau[13] + ked1_v2lapl2[0]*(ked2_v2sigma2[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*ked2_vsigma[0]*mgga_v3tau3[2] + 2*ked2_vsigma[0]*mgga_v3sigmatau2[7] + mgga_v3sigma2tau[10]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked2_v2sigma2[0]*mgga_v3tau3[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4tau4[2] +
    2*ked2_vsigma[0]*mgga_v4sigmatau3[9] + mgga_v4sigma2tau2[15]) +
    2*ked1_vlapl[0]*(ked2_v2sigma2[0]*mgga_v3lapltau2[1] + ked2_vsigma[0]*ked2_vsigma[0]*mgga_v4lapltau3[2] +
    2*ked2_vsigma[0]*mgga_v4sigmalapltau2[13] + mgga_v4sigma2lapltau[20]) + mgga_v4sigma2lapl2[15];
  v4sigma2lapl2[16] = ked2_v3sigma2lapl[0]*mgga_v2lapltau[1] + ked2_v2sigma2[0]*(ked2_vlapl[0]*mgga_v3lapltau2[2] +
    mgga_v3lapl2tau[3]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4lapltau3[3] + mgga_v4lapl2tau2[5]) +
    2*ked2_v2sigmalapl[0]*mgga_v3sigmalapltau[9] + 2*ked2_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3lapltau2[2] +
    ked2_vlapl[0]*mgga_v4sigmalapltau2[14] + mgga_v4sigmalapl2tau[15]) +
    ked1_vlapl[0]*(ked2_v3sigma2lapl[0]*mgga_v2tau2[1] + ked2_v2sigma2[0]*(ked2_vlapl[0]*mgga_v3tau3[2] +
    mgga_v3lapltau2[4]) + ked2_vsigma[0]*ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[3] + mgga_v4lapltau3[6]) +
    2*ked2_v2sigmalapl[0]*mgga_v3sigmatau2[7] + 2*ked2_vsigma[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] +
    ked2_vlapl[0]*mgga_v4sigmatau3[10] + mgga_v4sigmalapltau2[16]) + ked2_vlapl[0]*mgga_v4sigma2tau2[16] +
    mgga_v4sigma2lapltau[22]) + ked2_vlapl[0]*mgga_v4sigma2lapltau[21] + mgga_v4sigma2lapl2[16];
  v4sigma2lapl2[17] = ked2_v4sigma2lapl2[0]*mgga_vtau[1] + 2*ked2_v2sigmalapl[0]*ked2_v2sigmalapl[0]*mgga_v2tau2[2] +
    ked2_v2lapl2[0]*ked2_v2sigma2[0]*mgga_v2tau2[2] + 2*ked2_vlapl[0]*ked2_v3sigma2lapl[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_v2sigma2[0]*mgga_v3tau3[3] + 2*ked2_v3sigma2lapl[0]*mgga_v2lapltau[3] +
    2*ked2_vlapl[0]*ked2_v2sigma2[0]*mgga_v3lapltau2[5] + ked2_v2sigma2[0]*mgga_v3lapl2tau[5] +
    ked2_vsigma[0]*ked2_vsigma[0]*(ked2_v2lapl2[0]*mgga_v3tau3[3] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4tau4[4] +
    2*ked2_vlapl[0]*mgga_v4lapltau3[7] + mgga_v4lapl2tau2[8]) + 2*ked2_v3sigmalapl2[0]*mgga_v2sigmatau[5] +
    4*ked2_v2sigmalapl[0]*(ked2_vsigma[0]*mgga_v3lapltau2[5] + ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[3] +
    mgga_v3sigmatau2[8]) + mgga_v3sigmalapltau[11]) + 2*ked2_vsigma[0]*(ked2_v3sigmalapl2[0]*mgga_v2tau2[2] +
    ked2_v2lapl2[0]*mgga_v3sigmatau2[8] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmatau3[11] +
    2*ked2_vlapl[0]*mgga_v4sigmalapltau2[17] + mgga_v4sigmalapl2tau[17]) + ked2_v2lapl2[0]*mgga_v3sigma2tau[11] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigma2tau2[17] + 2*ked2_vlapl[0]*mgga_v4sigma2lapltau[23] +
    mgga_v4sigma2lapl2[17];
  v4sigma2lapltau[1] = 0;
  v4sigma2lapltau[2] = 0;
  v4sigma2lapltau[3] = 0;
  v4sigma2lapltau[4] = 0;
  v4sigma2lapltau[5] = 0;
  v4sigma2lapltau[6] = 0;
  v4sigma2lapltau[7] = 0;
  v4sigma2lapltau[8] = 0;
  v4sigma2lapltau[9] = 0;
  v4sigma2lapltau[10] = 0;
  v4sigma2lapltau[11] = 0;
  v4sigma2lapltau[12] = 0;
  v4sigma2lapltau[13] = 0;
  v4sigma2lapltau[14] = 0;
  v4sigma2lapltau[15] = 0;
  v4sigma2lapltau[16] = 0;
  v4sigma2lapltau[17] = 0;
  v4sigma2lapltau[18] = 0;
  v4sigma2lapltau[19] = 0;
  v4sigma2lapltau[20] = 0;
  v4sigma2lapltau[21] = 0;
  v4sigma2lapltau[22] = 0;
  v4sigma2lapltau[23] = 0;
  v4sigma2tau2[1] = 0;
  v4sigma2tau2[2] = 0;
  v4sigma2tau2[3] = 0;
  v4sigma2tau2[4] = 0;
  v4sigma2tau2[5] = 0;
  v4sigma2tau2[6] = 0;
  v4sigma2tau2[7] = 0;
  v4sigma2tau2[8] = 0;
  v4sigma2tau2[9] = 0;
  v4sigma2tau2[10] = 0;
  v4sigma2tau2[11] = 0;
  v4sigma2tau2[12] = 0;
  v4sigma2tau2[13] = 0;
  v4sigma2tau2[14] = 0;
  v4sigma2tau2[15] = 0;
  v4sigma2tau2[16] = 0;
  v4sigma2tau2[17] = 0;
  v4sigmalapl3[1] = ked1_v3sigmalapl2[0]*mgga_v2lapltau[2] + 2*ked1_vlapl[0]*ked1_v2sigmalapl[0]*mgga_v3lapltau2[3] +
    ked1_vlapl[0]*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4lapltau3[4] + 2*ked1_v2sigmalapl[0]*mgga_v3lapl2tau[2] +
    2*ked1_vlapl[0]*ked1_vsigma[0]*mgga_v4lapl2tau2[3] + ked1_vsigma[0]*mgga_v4lapl3tau[2] +
    ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3lapltau2[3] + mgga_v3sigmalapltau[2]) +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmalapltau2[3] + 2*ked1_vlapl[0]*mgga_v4sigmalapl2tau[2] +
    ked2_vlapl[0]*(ked1_v3sigmalapl2[0]*mgga_v2tau2[1] + 2*ked1_v2sigmalapl[0]*mgga_v3lapltau2[1] +
    ked1_vsigma[0]*mgga_v4lapl2tau2[1] + ked1_v2lapl2[0]*(ked1_vsigma[0]*mgga_v3tau3[1] + mgga_v3sigmatau2[1]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[1]) +
    2*ked1_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[1] + ked1_vsigma[0]*mgga_v4lapltau3[1] + mgga_v4sigmalapltau2[1])
    + mgga_v4sigmalapl2tau[1]) + mgga_v4sigmalapl3[1];
  v4sigmalapl3[2] = ked1_v2sigmalapl[0]*mgga_v3lapl2tau[4] + ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4lapl2tau2[6] +
    mgga_v4lapl3tau[4]) + ked1_vlapl[0]*mgga_v4sigmalapl2tau[4] + ked2_v2lapl2[0]*(ked1_v2sigmalapl[0]*mgga_v2tau2[1] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v3tau3[1] + mgga_v3lapltau2[1]) + ked1_vlapl[0]*mgga_v3sigmatau2[1] +
    mgga_v3sigmalapltau[1]) + ked2_vlapl[0]*ked2_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3tau3[2] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[2]) + ked1_vlapl[0]*mgga_v4sigmatau3[2] +
    mgga_v4sigmalapltau2[2]) + 2*ked2_vlapl[0]*(ked1_v2sigmalapl[0]*mgga_v3lapltau2[4] +
    ked1_vsigma[0]*(ked1_vlapl[0]*mgga_v4lapltau3[5] + mgga_v4lapl2tau2[4]) + ked1_vlapl[0]*mgga_v4sigmalapltau2[4] +
    mgga_v4sigmalapl2tau[3]) + mgga_v4sigmalapl3[2];
  v4sigmalapl3[3] = ked1_vsigma[0]*(3*ked2_v2lapl2[0]*mgga_v3lapltau2[4] + mgga_v4lapl3tau[6]) +
    ked2_v3lapl3[0]*(ked1_vsigma[0]*mgga_v2tau2[1] + mgga_v2sigmatau[1]) +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[3]) +
    3*ked2_v2lapl2[0]*mgga_v3sigmalapltau[3] + 3*ked2_vlapl[0]*ked2_vlapl[0]*(ked1_vsigma[0]*mgga_v4lapltau3[6] +
    mgga_v4sigmalapltau2[5]) + 3*ked2_vlapl[0]*(ked1_vsigma[0]*(ked2_v2lapl2[0]*mgga_v3tau3[2] + mgga_v4lapl2tau2[7]) +
    ked2_v2lapl2[0]*mgga_v3sigmatau2[2] + mgga_v4sigmalapl2tau[5]) + mgga_v4sigmalapl3[3];
  v4sigmalapl3[4] = ked1_v3lapl3[0]*mgga_v2sigmatau[2] + ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmatau3[4]
    + 3*ked1_v2lapl2[0]*mgga_v3sigmalapltau[4] + 3*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmalapltau2[6] +
    3*ked1_vlapl[0]*(ked1_v2lapl2[0]*mgga_v3sigmatau2[3] + mgga_v4sigmalapl2tau[6]) + mgga_v4sigmalapl3[4];
  v4sigmalapl3[5] = ked1_v2lapl2[0]*mgga_v3sigmalapltau[6] + ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmalapltau2[9] +
    2*ked1_vlapl[0]*mgga_v4sigmalapl2tau[8] + ked2_vlapl[0]*(ked1_v2lapl2[0]*mgga_v3sigmatau2[4] +
    ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4sigmatau3[5] + 2*ked1_vlapl[0]*mgga_v4sigmalapltau2[7] +
    mgga_v4sigmalapl2tau[7]) + mgga_v4sigmalapl3[5];
  v4sigmalapl3[6] = ked1_vlapl[0]*(ked2_v2lapl2[0]*mgga_v3sigmatau2[4] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmatau3[6] + 2*ked2_vlapl[0]*mgga_v4sigmalapltau2[10] +
    mgga_v4sigmalapl2tau[10]) + ked2_v2lapl2[0]*mgga_v3sigmalapltau[5] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmalapltau2[8] + 2*ked2_vlapl[0]*mgga_v4sigmalapl2tau[9] +
    mgga_v4sigmalapl3[6];
  v4sigmalapl3[7] = ked2_v3lapl3[0]*mgga_v2sigmatau[3] + ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmatau3[7]
    + 3*ked2_v2lapl2[0]*mgga_v3sigmalapltau[7] + 3*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmalapltau2[11] +
    3*ked2_vlapl[0]*(ked2_v2lapl2[0]*mgga_v3sigmatau2[5] + mgga_v4sigmalapl2tau[11]) + mgga_v4sigmalapl3[7];
  v4sigmalapl3[8] = ked2_vsigma[0]*(3*ked1_v2lapl2[0]*mgga_v3lapltau2[1] + mgga_v4lapl3tau[1]) +
    ked1_v3lapl3[0]*(ked2_vsigma[0]*mgga_v2tau2[1] + mgga_v2sigmatau[4]) +
    ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[1] + mgga_v4sigmatau3[8]) +
    3*ked1_v2lapl2[0]*mgga_v3sigmalapltau[8] + 3*ked1_vlapl[0]*ked1_vlapl[0]*(ked2_vsigma[0]*mgga_v4lapltau3[1] +
    mgga_v4sigmalapltau2[12]) + 3*ked1_vlapl[0]*(ked2_vsigma[0]*(ked1_v2lapl2[0]*mgga_v3tau3[1] + mgga_v4lapl2tau2[1])
    + ked1_v2lapl2[0]*mgga_v3sigmatau2[6] + mgga_v4sigmalapl2tau[12]) + mgga_v4sigmalapl3[8];
  v4sigmalapl3[9] = ked2_v2sigmalapl[0]*mgga_v3lapl2tau[1] + ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4lapl2tau2[2] +
    mgga_v4lapl3tau[3]) + ked1_v2lapl2[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[1] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v3tau3[2] + mgga_v3lapltau2[4]) + ked2_vlapl[0]*mgga_v3sigmatau2[7] +
    mgga_v3sigmalapltau[10]) + ked1_vlapl[0]*ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[1] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4tau4[2] + mgga_v4lapltau3[5]) + ked2_vlapl[0]*mgga_v4sigmatau3[9] +
    mgga_v4sigmalapltau2[15]) + 2*ked1_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3lapltau2[1] +
    ked2_vsigma[0]*(ked2_vlapl[0]*mgga_v4lapltau3[2] + mgga_v4lapl2tau2[4]) + ked2_vlapl[0]*mgga_v4sigmalapltau2[13] +
    mgga_v4sigmalapl2tau[14]) + ked2_vlapl[0]*mgga_v4sigmalapl2tau[13] + mgga_v4sigmalapl3[9];
  v4sigmalapl3[10] = ked2_v3sigmalapl2[0]*mgga_v2lapltau[1] + 2*ked2_vlapl[0]*ked2_v2sigmalapl[0]*mgga_v3lapltau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4lapltau3[3] + 2*ked2_v2sigmalapl[0]*mgga_v3lapl2tau[3] +
    2*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4lapl2tau2[5] + ked2_vsigma[0]*mgga_v4lapl3tau[5] +
    ked1_vlapl[0]*(ked2_v3sigmalapl2[0]*mgga_v2tau2[1] + 2*ked2_v2sigmalapl[0]*mgga_v3lapltau2[4] +
    ked2_vsigma[0]*mgga_v4lapl2tau2[7] + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3tau3[2] + mgga_v3sigmatau2[7]) +
    ked2_vlapl[0]*ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v4tau4[3] + mgga_v4sigmatau3[10]) +
    2*ked2_vlapl[0]*(ked2_v2sigmalapl[0]*mgga_v3tau3[2] + ked2_vsigma[0]*mgga_v4lapltau3[6] + mgga_v4sigmalapltau2[16])
    + mgga_v4sigmalapl2tau[16]) + ked2_v2lapl2[0]*(ked2_vsigma[0]*mgga_v3lapltau2[2] + mgga_v3sigmalapltau[9]) +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmalapltau2[14] + 2*ked2_vlapl[0]*mgga_v4sigmalapl2tau[15] +
    mgga_v4sigmalapl3[10];
  v4sigmalapl3[11] = ked2_v4sigmalapl3[0]*mgga_vtau[1] + 3*ked2_vlapl[0]*ked2_v3sigmalapl2[0]*mgga_v2tau2[2] +
    3*ked2_vlapl[0]*ked2_vlapl[0]*ked2_v2sigmalapl[0]*mgga_v3tau3[3] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4tau4[4] + 3*ked2_v3sigmalapl2[0]*mgga_v2lapltau[3]
    + 6*ked2_vlapl[0]*ked2_v2sigmalapl[0]*mgga_v3lapltau2[5] +
    3*ked2_vlapl[0]*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4lapltau3[7] + 3*ked2_v2sigmalapl[0]*mgga_v3lapl2tau[5] +
    3*ked2_vlapl[0]*ked2_vsigma[0]*mgga_v4lapl2tau2[8] + ked2_vsigma[0]*mgga_v4lapl3tau[7] +
    ked2_v3lapl3[0]*(ked2_vsigma[0]*mgga_v2tau2[2] + mgga_v2sigmatau[5]) +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmatau3[11] +
    3*ked2_v2lapl2[0]*(ked2_v2sigmalapl[0]*mgga_v2tau2[2] + ked2_vsigma[0]*mgga_v3lapltau2[5] +
    ked2_vlapl[0]*(ked2_vsigma[0]*mgga_v3tau3[3] + mgga_v3sigmatau2[8]) + mgga_v3sigmalapltau[11]) +
    3*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4sigmalapltau2[17] + 3*ked2_vlapl[0]*mgga_v4sigmalapl2tau[17] +
    mgga_v4sigmalapl3[11];
  v4sigmalapl2tau[1] = 0;
  v4sigmalapl2tau[2] = 0;
  v4sigmalapl2tau[3] = 0;
  v4sigmalapl2tau[4] = 0;
  v4sigmalapl2tau[5] = 0;
  v4sigmalapl2tau[6] = 0;
  v4sigmalapl2tau[7] = 0;
  v4sigmalapl2tau[8] = 0;
  v4sigmalapl2tau[9] = 0;
  v4sigmalapl2tau[10] = 0;
  v4sigmalapl2tau[11] = 0;
  v4sigmalapl2tau[12] = 0;
  v4sigmalapl2tau[13] = 0;
  v4sigmalapl2tau[14] = 0;
  v4sigmalapl2tau[15] = 0;
  v4sigmalapl2tau[16] = 0;
  v4sigmalapl2tau[17] = 0;
  v4sigmalapltau2[1] = 0;
  v4sigmalapltau2[2] = 0;
  v4sigmalapltau2[3] = 0;
  v4sigmalapltau2[4] = 0;
  v4sigmalapltau2[5] = 0;
  v4sigmalapltau2[6] = 0;
  v4sigmalapltau2[7] = 0;
  v4sigmalapltau2[8] = 0;
  v4sigmalapltau2[9] = 0;
  v4sigmalapltau2[10] = 0;
  v4sigmalapltau2[11] = 0;
  v4sigmalapltau2[12] = 0;
  v4sigmalapltau2[13] = 0;
  v4sigmalapltau2[14] = 0;
  v4sigmalapltau2[15] = 0;
  v4sigmalapltau2[16] = 0;
  v4sigmalapltau2[17] = 0;
  v4sigmatau3[1] = 0;
  v4sigmatau3[2] = 0;
  v4sigmatau3[3] = 0;
  v4sigmatau3[4] = 0;
  v4sigmatau3[5] = 0;
  v4sigmatau3[6] = 0;
  v4sigmatau3[7] = 0;
  v4sigmatau3[8] = 0;
  v4sigmatau3[9] = 0;
  v4sigmatau3[10] = 0;
  v4sigmatau3[11] = 0;
  v4lapl4[1] = ked1_v3lapl3[0]*mgga_v2lapltau[2] + ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4lapltau3[4] +
    3*ked1_v2lapl2[0]*mgga_v3lapl2tau[2] + 3*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4lapl2tau2[3] +
    3*ked1_vlapl[0]*(ked1_v2lapl2[0]*mgga_v3lapltau2[3] + mgga_v4lapl3tau[2]) +
    ked2_vlapl[0]*(ked1_v3lapl3[0]*mgga_v2tau2[1] + ked1_vlapl[0]*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4tau4[1] +
    3*ked1_v2lapl2[0]*mgga_v3lapltau2[1] + 3*ked1_vlapl[0]*ked1_vlapl[0]*mgga_v4lapltau3[1] +
    3*ked1_vlapl[0]*(ked1_v2lapl2[0]*mgga_v3tau3[1] + mgga_v4lapl2tau2[1]) + mgga_v4lapl3tau[1]) + mgga_v4lapl4[1];
  v4lapl4[2] = ked1_v2lapl2[0]*(ked2_v2lapl2[0]*mgga_v2tau2[1] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3tau3[2] +
    2*ked2_vlapl[0]*mgga_v3lapltau2[4] + mgga_v3lapl2tau[4]) +
    ked1_vlapl[0]*ked1_vlapl[0]*(ked2_v2lapl2[0]*mgga_v3tau3[1] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4tau4[2] +
    2*ked2_vlapl[0]*mgga_v4lapltau3[5] + mgga_v4lapl2tau2[6]) + 2*ked1_vlapl[0]*(ked2_v2lapl2[0]*mgga_v3lapltau2[1] +
    ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapltau3[2] + 2*ked2_vlapl[0]*mgga_v4lapl2tau2[4] + mgga_v4lapl3tau[4]) +
    ked2_v2lapl2[0]*mgga_v3lapl2tau[1] + ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapl2tau2[2] +
    2*ked2_vlapl[0]*mgga_v4lapl3tau[3] + mgga_v4lapl4[2];
  v4lapl4[3] = ked1_vlapl[0]*(ked2_v3lapl3[0]*mgga_v2tau2[1] + ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4tau4[3]
    + 3*ked2_v2lapl2[0]*mgga_v3lapltau2[4] + 3*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapltau3[6] +
    3*ked2_vlapl[0]*(ked2_v2lapl2[0]*mgga_v3tau3[2] + mgga_v4lapl2tau2[7]) + mgga_v4lapl3tau[6]) +
    ked2_v3lapl3[0]*mgga_v2lapltau[1] + ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapltau3[3] +
    3*ked2_v2lapl2[0]*mgga_v3lapl2tau[3] + 3*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapl2tau2[5] +
    3*ked2_vlapl[0]*(ked2_v2lapl2[0]*mgga_v3lapltau2[2] + mgga_v4lapl3tau[5]) + mgga_v4lapl4[3];
  v4lapl4[4] = ked2_v4lapl4[0]*mgga_vtau[1] + 3*ked2_v2lapl2[0]*ked2_v2lapl2[0]*mgga_v2tau2[2] +
    ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4tau4[4] + 4*ked2_v3lapl3[0]*mgga_v2lapltau[3] +
    4*ked2_vlapl[0]*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapltau3[7] +
    6*ked2_v2lapl2[0]*(ked2_vlapl[0]*ked2_vlapl[0]*mgga_v3tau3[3] + 2*ked2_vlapl[0]*mgga_v3lapltau2[5] +
    mgga_v3lapl2tau[5]) + 6*ked2_vlapl[0]*ked2_vlapl[0]*mgga_v4lapl2tau2[8] +
    4*ked2_vlapl[0]*(ked2_v3lapl3[0]*mgga_v2tau2[2] + mgga_v4lapl3tau[7]) + mgga_v4lapl4[4];
  v4lapl3tau[1] = 0;
  v4lapl3tau[2] = 0;
  v4lapl3tau[3] = 0;
  v4lapl3tau[4] = 0;
  v4lapl3tau[5] = 0;
  v4lapl3tau[6] = 0;
  v4lapl3tau[7] = 0;
  v4lapl2tau2[1] = 0;
  v4lapl2tau2[2] = 0;
  v4lapl2tau2[3] = 0;
  v4lapl2tau2[4] = 0;
  v4lapl2tau2[5] = 0;
  v4lapl2tau2[6] = 0;
  v4lapl2tau2[7] = 0;
  v4lapl2tau2[8] = 0;
  v4lapltau3[1] = 0;
  v4lapltau3[2] = 0;
  v4lapltau3[3] = 0;
  v4lapltau3[4] = 0;
  v4lapltau3[5] = 0;
  v4lapltau3[6] = 0;
  v4lapltau3[7] = 0;
  v4tau4[1] = 0;
  v4tau4[2] = 0;
  v4tau4[3] = 0;
  v4tau4[4] = 0;
}
