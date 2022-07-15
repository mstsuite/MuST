(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define mgga_x_revtpss_params
$include "mgga_x_tpss.mpl"

reg_c := 3:
reg_d := 1.475:

reg_f_a :=  a -> (1 - a)^3/(1 + (reg_d*a)^2)^(3/2):

(* Eq. (12). Note that alpha = 0 => t = x^2/8 *)
reg_f := (x, u, t) ->
  tpss_f(x, u, t) + reg_f_a(tpss_alpha(x, t))*exp(-reg_c*X2S^2*x^2)*(tpss_f(x, u, x^2/8) - tpss_f(x, u, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(reg_f, rs, z, xs0, xs1, u0, u1, t0, t1):
