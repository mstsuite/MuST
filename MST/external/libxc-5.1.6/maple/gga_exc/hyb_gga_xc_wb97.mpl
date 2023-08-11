(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_xc_wb97_params *params;

  assert(p->params != NULL);
  params = (gga_xc_wb97_params * )(p->params);
*)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$include "lda_x_erf.mpl"

$include "b97.mpl"

wb97_x := (rs, z, xs0, xs1) ->
  my_piecewise3(screen_dens_zeta(rs,  z), 0, (1+z)/2 * lda_x_erf_spin(rs*(2/(1 + z))^(1/3),  1) * b97_g(0.004, params_a_c_x, xs0))
+ my_piecewise3(screen_dens_zeta(rs, -z), 0, (1-z)/2 * lda_x_erf_spin(rs*(2/(1 - z))^(1/3),  1) * b97_g(0.004, params_a_c_x, xs1)):

f := (rs, z, xt, xs0, xs1) ->
  wb97_x(rs, z, xs0, xs1) +
  b97_f(f_pw, 0.2, params_a_c_ss, 0.006, params_a_c_ab, rs, z, xs0, xs1):

