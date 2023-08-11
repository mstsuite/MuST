(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_ft97_params *params;

  assert(p->params != NULL);
  params = (gga_x_ft97_params * )(p->params);
*)

ft97_beta := (rs, z, xs) -> params_a_beta0
  + params_a_beta1*sigma_spin(rs, z, xs)/(params_a_beta2 + sigma_spin(rs, z, xs)):

ft97_fx := (rs, z, xs) -> 1 + ft97_beta(rs, z, xs)*xs^2 /
  (X_FACTOR_C*sqrt(1 + 9*xs^2*ft97_beta(rs, z, xs)^2*arcsinh(xs^2)^2)):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange_nsp(ft97_fx, rs, zeta, xs0, xs1):
