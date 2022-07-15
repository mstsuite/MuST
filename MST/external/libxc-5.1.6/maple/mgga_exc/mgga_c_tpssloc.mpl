(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "gga_c_pbeloc.mpl"

params_a_C0_c := [0.35, 0.87, 0.50, 2.26]:
params_a_d    := 4.5:
$include "tpss_c.mpl"

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  + tpss_f(f_pbe, rs, z, xt, xs0, xs1, ts0, ts1):
