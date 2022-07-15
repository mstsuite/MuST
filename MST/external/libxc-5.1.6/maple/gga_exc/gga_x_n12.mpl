(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_n12_params *params;

  assert(p->params != NULL);
  params = (gga_x_n12_params * )(p->params);
*)

n12_omega_x := 2.5:
n12_gamma_x := 0.004:

n12_rss := (rs, z) -> rs * 2^(1/3) * opz_pow_n(z,-1/3):

n12_vx := rs -> 1/(1 + (1/(RS_FACTOR*n12_omega_x))*rs):
n12_ux := x -> n12_gamma_x*x^2/(1 + n12_gamma_x*x^2):

n12_FN12 := (rs, z, x) ->
  + add(params_a_CC_0_[i+1]*n12_ux(x)^i, i=0..3)
  + add(params_a_CC_1_[i+1]*n12_ux(x)^i, i=0..3) * n12_vx(n12_rss(rs, z))
  + add(params_a_CC_2_[i+1]*n12_ux(x)^i, i=0..3) * n12_vx(n12_rss(rs, z))^2
  + add(params_a_CC_3_[i+1]*n12_ux(x)^i, i=0..3) * n12_vx(n12_rss(rs, z))^3:

f  := (rs, z, xt, xs0, xs1) -> gga_exchange_nsp(n12_FN12, rs, z, xs0, xs1):
