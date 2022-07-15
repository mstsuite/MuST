(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

gea2_s := x -> X2S*x:
gea2_q := u -> X2S^2*u:

gea2_f0 := (s, q) -> 1 + 5/27*s^2 + 20/9*q:
gea2_f := (x, u) ->
  gea2_f0(gea2_s(x), gea2_q(u)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(gea2_f, rs, z, xs0, xs1, u0, u1):
