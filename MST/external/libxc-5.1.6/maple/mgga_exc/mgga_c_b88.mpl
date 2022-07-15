(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define gga_x_b86_mgc_params
$include "gga_x_b86.mpl"

cab := 0.63:
css := 0.96:

(* Equation 50, same-spin correlation *)
b88_css := (rs, z, xs, ts) ->
  - 0.01 * (1 + z)^(8/3) * 2^(-8/3) * n_total(rs)^(5/3) * (2*ts - xs^2/4)
  * b88_zss(css, b86_f, rs, z, xs)^4 * (
    1 - 2*log(1 + b88_zss(css, b86_f, rs, z, xs)/2)
      / b88_zss(css, b86_f, rs, z, xs)
  ):

(* Same-spin correlation overall *)
b88_par := (rs, z, xs0, xs1, ts0, ts1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, b88_css(rs, z_thr( z), xs0, ts0))
  + my_piecewise3(screen_dens(rs, -z), 0, b88_css(rs, z_thr(-z), xs1, ts1)):

(* Equation 49, opposite-spin correlation *)
b88_cab := (rs, z, xs0, xs1) ->
  - 0.8 * (1 - z^2)/4 * n_total(rs)
  * b88_zab(cab, b86_f, rs, z, xs0, xs1) * (
      b88_zab(cab, b86_f, rs, z, xs0, xs1) - log(1 + b88_zab(cab, b86_f, rs, z, xs0, xs1))
  ):

(* Whole functional *)
b88_c_f := (rs, z, xs0, xs1, ts0, ts1) ->
  + b88_cab(rs,  z, xs0, xs1)
  + b88_par(rs,  z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  b88_c_f(rs, z, xs0, xs1, ts0, ts1):
