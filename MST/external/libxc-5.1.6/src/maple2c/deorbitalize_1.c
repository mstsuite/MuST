vrho[0] = ked1_vrho[0]*mgga_vtau[0] + mgga_vrho[0];
vsigma[0] = ked1_vsigma[0]*mgga_vtau[0] + mgga_vsigma[0];
vlapl[0] = ked1_vlapl[0]*mgga_vtau[0] + mgga_vlapl[0];
vtau[0] = 0;

if(func->nspin == XC_POLARIZED){
  vrho[1] = ked2_vrho[0]*mgga_vtau[1] + mgga_vrho[1];
  vsigma[1] = mgga_vsigma[1];
  vsigma[2] = ked2_vsigma[0]*mgga_vtau[1] + mgga_vsigma[2];
  vlapl[1] = ked2_vlapl[0]*mgga_vtau[1] + mgga_vlapl[1];
  vtau[1] = 0;
}
