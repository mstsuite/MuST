(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_pgslb_params *params;

  assert(p->params != NULL);
  params = (mgga_k_pgslb_params * )(p->params);
*)

# Equation (4) and (8)
pgslb_f0 := (s, q) -> 5/3*s^2 + exp(-params_a_pgslb_mu * s^2) + params_a_pgslb_beta*q^2:
pgslb_f := (x, u) -> pgslb_f0(X2S*x, X2S^2*u):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(pgslb_f, rs, z, xs0, xs1, u0, u1):

