(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* replace: "br89_x\(" -> "xc_mgga_x_br89_get_x(" *)

$include "mgga_x_br89.mpl"

mbrxh_a1 := 0.23432:
mbrxh_a2 := 0.089:
mbrxh_a3 := 0.0053:
params_a_at := 0:

(* new definition of Q. The rest of the functional remains the same *)
br89_Q := (x, u, t) ->
      mbrxh_a1*(2*t) - K_FACTOR_C + mbrxh_a2*x^2 + mbrxh_a3*x^4: