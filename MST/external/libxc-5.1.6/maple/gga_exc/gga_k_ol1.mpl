(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

ol1_f  := x -> 1 + (x^2/72 + 0.00677*2^(1/3)*x)/K_FACTOR_C:

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(ol1_f, rs, zeta, xs0, xs1):
