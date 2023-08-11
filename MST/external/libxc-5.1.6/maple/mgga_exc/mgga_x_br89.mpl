(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_br89_params *params;

  assert(p->params != NULL);
  params = (mgga_x_br89_params * ) (p->params);
*)

(* replace: "br89_x\(" -> "xc_mgga_x_br89_get_x(" *)

(* This is the derivative of f = x*exp(-2.0/3.0*x)/(x - 2) = y = 2*Pi^(2/3)/(3*Q) *)
br89_aux_dfdx := x -> -2/3 * exp(-2*x/3) * (x^2 - 2*x + 3) / (x - 2)^2:

`diff/br89_x` := proc(Q, g)
  -2/3 * Pi^(2/3) * diff(Q, g)/(Q^2 * br89_aux_dfdx(br89_x(Q)))
end proc:

br89_Q := (x, u, t) ->
  (u - 4*params_a_gamma*t + params_a_gamma*x^2/2)/6:

br89_min_Q := 5.0e-13:
br89_cQ := Q -> my_piecewise3(abs(Q) < br89_min_Q,
  my_piecewise3(Q > 0, br89_min_Q, -br89_min_Q), Q):

br89_v := x ->
  -2*Pi^(1/3)/X_FACTOR_C * exp(x/3)*(1 - exp(-x)*(1 + x/2))/x:

br89_mx := Q -> br89_x(Q):

br89_f := (x, u, t) ->
  - br89_v(br89_mx(br89_cQ(br89_Q(x, u, t))))/2 *
  (1 + params_a_at*mgga_series_w([0, 1, 0, -2, 0, 1], 6, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(br89_f, rs, z, xs0, xs1, u0, u1, t0, t1):
