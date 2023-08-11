(*
 Copyright (C) 2019 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "hyb_mgga_x_pjs18.mpl"

(* This expression (10) has \tilde A, and not A *)
js18_f_SR := (rs, z, x, t) -> tm_w(x, t)*js18_DME_SR(rs, z, x, t)
  + (1 - tm_w(x, t))*attenuation_erf(a_cnst*rs/opz_pow_n(z,1/3))*tm_fx_SC(x, t):

js18_f := (rs, z, x, u, t) -> -p_a_hyb_coeff_0_*js18_f_SR(rs, z, x, t) + tm_f(x, u, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange_nsp(js18_f, rs, z, xs0, xs1, u0, u1, t0, t1):
