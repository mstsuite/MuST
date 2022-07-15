(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

gamma_ss := 0.2:
cc_ss    := [0.0136823, 0.268920, -0.550769,  1.03947, 0.0]:

gamma_ab := 0.006:
cc_ab    := [0.836897,  1.72051,  -2.78498,  -4.57504, 0.0]:

$include "lda_c_vwn.mpl"
$include "b97.mpl"

f := (rs, z, xt, xs0, xs1) ->
  b97_f(f_vwn, gamma_ss, cc_ss, gamma_ab, cc_ab,
        rs, z, xs0, xs1):