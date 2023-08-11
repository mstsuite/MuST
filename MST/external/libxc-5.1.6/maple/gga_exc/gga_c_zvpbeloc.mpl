(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "gga_c_zvpbeint.mpl"
$include "gga_c_pbeloc.mpl"

params_a_alpha := 0.5:
params_a_omega := 2:

(* redefine nu of zbpbeint *)
zvpbeint_nu := (rs, z, t) ->
  2*(4/(3*Pi^2))^(1/18) * rs^(1/3):

(* Note that f_pbe here is, in fact, pbeloc *)
f  := (rs, z, xt, xs0, xs1) ->
  zvpbeint_ff(rs, z, 0) * f_pbe(rs, z, xt, xs0, xs1):
