(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

params_a_n := 19:

params_a_a := [
   13/12, 7/6,  8/6,  9/6, 10/6, 17/12, 9/6, 10/6,
   11/6, 10/6, 11/6, 12/6, 10/6, 11/6, 12/6,  7/6,
    8/6,  9/6, 10/6.0
]:

params_a_b := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]:
params_a_c := [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0]:
params_a_d := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]:

params_a_omega := [
   +0.678831e+00, -0.175821e+01, +0.127676e+01, -0.160789e+01, +0.365610e+00, -0.181327e+00,
   +0.146973e+00, +0.147141e+00, -0.716917e-01, -0.407167e-01, +0.214625e-01, -0.768156e-03,
   +0.310377e-01, -0.720326e-01, +0.446562e-01, -0.266802e+00, +0.150822e+01, -0.194515e+01,
   +0.679078e+00
]:

$include "th.mpl"

f := (rs, z, xt, xs0, xs1) -> f_th(rs, z, xt, xs0, xs1):
