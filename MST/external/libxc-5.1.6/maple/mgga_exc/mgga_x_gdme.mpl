(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_gdme_params *params;

  assert(p->params != NULL);
  params = (mgga_x_gdme_params * )(p->params);
*)

gdme_at := (params_a_AA + 3/5*params_a_BB)*2^(1/3)/(X_FACTOR_C*(3*Pi^2)^(2/3)):
gdme_bt := params_a_BB/(X_FACTOR_C*2^(1/3)*(3*Pi^2)^(4/3)):

gdme_f := (x, u, t) -> gdme_at + gdme_bt*((params_a_a^2 - params_a_a + 1/2)*u - 2*t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(gdme_f, rs, z, xs0, xs1, u0, u1, t0, t1):
