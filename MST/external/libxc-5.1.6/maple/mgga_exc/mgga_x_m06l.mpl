(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_m06l_params *params;

  assert(p->params != NULL);
  params = (mgga_x_m06l_params * )(p->params);
*)

$define gga_x_pbe_params
$include "gga_x_pbe.mpl"
$include "gvt4.mpl"

m06_alpha  := 0.00186726:
m06_coeff_d := params_a_d:

(* there is a factor if 2 in the definition of z, as in Theor. Chem. Account 120, 215 (2008) *)
(* A MINUS was missing in Eq. (7) of the paper *)

m06_f := (x, u, t) ->
  + pbe_f(x)*mgga_series_w(params_a_a, 12, t)
  + gtv4(m06_alpha, m06_coeff_d, x, 2*(t - K_FACTOR_C)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(m06_f, rs, z, xs0, xs1, u0, u1, t0, t1):
