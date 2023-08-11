(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_m11_params *params;

  assert(p->params != NULL);
  params = (mgga_x_m11_params * ) (p->params);
*)

$include "mgga_x_m08.mpl"
$include "lda_x_erf.mpl"

m11_f := (rs, z, x, u, t) ->
   attenuation_erf(a_cnst*rs/opz_pow_n(z,1/3)) * m08_f0(params_a_a, params_a_b, x, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange_nsp(m11_f, rs, z, xs0, xs1, u0, u1, t0, t1):
