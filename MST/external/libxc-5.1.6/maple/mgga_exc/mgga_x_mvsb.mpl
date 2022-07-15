(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_mvsb_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mvsb_params * ) (p->params);
*)

$include "mgga_x_mvs.mpl"

mvsb_beta := (t, x) -> mvs_alpha(t, x)*K_FACTOR_C/(t - K_FACTOR_C):

mvsb_f := (x, u, t) -> (1 + params_a_k0*mvs_fa(mvsb_beta(t,x)))
       / (1 + params_a_b*(X2S*x)^4)^(1/8):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(mvsb_f, rs, z, xs0, xs1, u0, u1, t0, t1):
