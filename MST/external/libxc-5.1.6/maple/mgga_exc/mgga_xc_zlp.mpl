(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

zlp_c := 0.828432   *RS_FACTOR:
zlp_d := 2.15509e-2 *RS_FACTOR:
zlp_k := 2.047107e-3*RS_FACTOR:

zlp_f := (rs, z, xt, us0, us1) ->
  - (zlp_c + zlp_d*t_vw(z, xt, us0, us1))
  * (1 - zlp_k*log(1 + rs/zlp_k)/rs)/rs:

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  zlp_f(rs, z, xt, us0, us1):
