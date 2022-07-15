(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_m08_params *params;

  assert(p->params != NULL);
  params = (mgga_x_m08_params * ) (p->params);
*)

params_a_rpbe_kappa := 0.552:
params_a_rpbe_mu    := MU_GE:
$include "gga_x_rpbe.mpl"

params_a_kappa := KAPPA_PBE:
params_a_mu    := 0.21951:
$include "gga_x_pbe.mpl"

m08_f0 := (a, b, x, t) ->
  + pbe_f(x) *mgga_series_w(a, 12, t)
  + rpbe_f(x)*mgga_series_w(b, 12, t):

m08_f := (x, u, t) -> m08_f0(params_a_a, params_a_b, x, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(m08_f, rs, z, xs0, xs1, u0, u1, t0, t1):
