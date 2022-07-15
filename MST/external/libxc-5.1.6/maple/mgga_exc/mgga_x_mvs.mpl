(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_mvs_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mvs_params * ) (p->params);
*)

(* equation 10 *)
mvs_fa := a -> (1 - a) / ((1 + params_a_e1*a^2)^2 + params_a_c1*a^4)^(1/4):

(* alpha *)
mvs_alpha := (t,x) -> (t - x^2/8)/K_FACTOR_C:

(* eq 7 *)
mvs_f := (x, u, t) -> (1 + params_a_k0*mvs_fa(mvs_alpha(t,x)))
      / (1 + params_a_b*(X2S*x)^4)^(1/8):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(mvs_f, rs, z, xs0, xs1, u0, u1, t0, t1):