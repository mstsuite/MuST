(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

gvt4_alpha  := 0.00186726:
gvt4_coeff_d := [-9.800683e-01, -3.556788e-03, 6.250326e-03, -2.354518e-05, -1.282732e-04, 3.574822e-04]:

$include "gvt4.mpl"

gvt4_f := (x, u, t) -> -gtv4(gvt4_alpha, gvt4_coeff_d, x, 2*(t - K_FACTOR_C))/X_FACTOR_C:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(gvt4_f, rs, z, xs0, xs1, u0, u1, t0, t1):
