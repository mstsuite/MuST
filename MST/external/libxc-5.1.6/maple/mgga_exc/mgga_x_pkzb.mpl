(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

pkzb_a2 := 146/2025:
pkzb_a3 := -73/405:
pkzb_a4 := 0.131957187845257783631757384393: (* DD + 100/(81*81*kappa) *)

pkzb_kappa := 0.804:

pkzb_xx := (p, qt) -> MU_GE*p + pkzb_a2*qt*qt + pkzb_a3*qt*p + pkzb_a4*p*p:

pkzb_f := (x, u, t) -> 1 + pkzb_kappa - \
  pkzb_kappa^2/(pkzb_kappa + pkzb_xx(X2S^2*x^2, 6*X2S^2*t - 9/20 - X2S^2*x^2/12)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(pkzb_f, rs, z, xs0, xs1, u0, u1, t0, t1):