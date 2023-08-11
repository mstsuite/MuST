(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "gga_x_pw86.mpl"

params_a_aa := 2.208:
params_a_bb := 9.27:
params_a_cc := 0.2:

f := (rs, z, xt, xs0, xs1) -> gga_kinetic(pw86_f, rs, z, xs0, xs1):
