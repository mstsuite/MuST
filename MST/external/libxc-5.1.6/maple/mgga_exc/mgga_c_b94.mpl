(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_c_b94_params *params;

  assert(p->params != NULL);
  params = (mgga_c_b94_params * ) (p->params);
*)

(* replace: "br89_x\(" -> "xc_mgga_x_br89_get_x(" *)

$define mgga_x_br89_params
$include "mgga_x_br89.mpl"

(* This is a fake parameter in libxc *)
params_a_at := 0:

(* Equation 9, same-spin correlation *)
b94_css := (rs, z, xs, us, ts) ->
  - 0.01 * (1 + z)^(8/3) * 2^(-8/3) * n_total(rs)^(5/3) * (2*ts - xs^2/4)
  * b94_zss(params_a_css, br89_f, rs, z, xs, us, ts)^4 * (
    1 - 2*log(1 + b94_zss(params_a_css, br89_f, rs, z, xs, us, ts)/2)
      / b94_zss(params_a_css, br89_f, rs, z, xs, us, ts)
  ):

(* Same-spin correlation overall *)
b94_par := (rs, z, xs0, xs1, us0, us1, ts0, ts1) ->
  + my_piecewise3(screen_dens(rs,  z), 0, b94_css(rs, z_thr( z), xs0, us0, ts0))
  + my_piecewise3(screen_dens(rs, -z), 0, b94_css(rs, z_thr(-z), xs1, us1, ts1)):

(* Equation 8, opposite-spin correlation *)
b94_cab := (rs, z, xs0, xs1, us0, us1, ts0, ts1) ->
  - 0.8 * (1 - z^2)/4 * n_total(rs)
  * b94_zab(params_a_cab, br89_f, rs, z, xs0, xs1, us0, us1, ts0, ts1) * (
     b94_zab(params_a_cab, br89_f, rs, z, xs0, xs1, us0, us1, ts0, ts1)
     - log(1 + b94_zab(params_a_cab, br89_f, rs, z, xs0, xs1, us0, us1, ts0, ts1))
  ):

(* Whole functional *)
b94_c_f := (rs, z, xs0, xs1, us0, us1, ts0, ts1) ->
  + b94_cab(rs,  z, xs0, xs1, us0, us1, ts0, ts1)
  + b94_par(rs,  z, xs0, xs1, us0, us1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  b94_c_f(rs, z, xs0, xs1, us0, us1, ts0, ts1):
