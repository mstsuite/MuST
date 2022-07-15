(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_lta_params *params;

  assert(p->params != NULL);
  params = (mgga_x_lta_params * )(p->params);
*)

lta_f := (x, u, t) -> (t/K_FACTOR_C)^(4*params_a_ltafrac/5):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(lta_f, rs, z, xs0, xs1, u0, u1, t0, t1):
