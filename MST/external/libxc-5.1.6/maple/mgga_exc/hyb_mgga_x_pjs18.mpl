(*
 Copyright (C) 2019 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_x_tm.mpl"
$include "lda_x_erf.mpl"

(* See text after eq 7  *)
js18_H := x -> (2*tm_lambda - 1)^2 * tm_p(x):
(* See text after eq 7 *)
js18_G := (x, t) ->
  (3*(tm_lambda^2 - tm_lambda + 1/2)*(t - K_FACTOR_C - x^2/72) - (t - K_FACTOR_C)
   + 7/18*(2*tm_lambda - 1)^2*x^2)/K_FACTOR_C:

(* not to run into problems with erfs and exp *)
js18_A := (rs, z, x) -> m_max(1e-10, a_cnst*rs/(tm_f0(x)*opz_pow_n(z,1/3))):

js18_DME_SR := (rs, z, x, t) ->
  + attenuation_erf   (js18_A(rs, z, x))/tm_f0(x)^2
  + attenuation_erf_f2(js18_A(rs, z, x))*7*js18_G(x, t)/(9*tm_f0(x)^4)
  + attenuation_erf_f3(js18_A(rs, z, x))*245*js18_H(x)/(54*tm_f0(x)^4):

pjs18_f := (rs, z, x, u, t) -> js18_DME_SR(rs, z, x, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange_nsp(pjs18_f, rs, z, xs0, xs1, u0, u1, t0, t1):
