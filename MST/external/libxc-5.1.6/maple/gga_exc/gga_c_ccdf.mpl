(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_ccdf_params *params;

  assert(p->params != NULL);
  params = (gga_c_ccdf_params * )(p->params);
*)

(* Equation (26) *)
f_ccdf  := (rs, z, xt, xs0, xs1) ->
  params_a_c1 / (1 + params_a_c2*n_total(rs)^(-1/3)) * (1 - params_a_c3 / (1 + exp(-params_a_c4*(2^(1/3)*X2S*xt - params_a_c5)))):

f  := (rs, z, xt, xs0, xs1) -> f_ccdf(rs, z, xt, xs0, xs1):
