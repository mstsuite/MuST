(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

herman_c1 := 0.003:

herman_f := x -> 1 + herman_c1/X_FACTOR_C * x^2:

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(herman_f, rs, zeta, xs0, xs1):
