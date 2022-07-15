(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

gea4_s := x -> X2S*x:
gea4_q := u -> X2S^2*u:

gea4_f0 := (s, q) -> 1 + 5/27*s^2 + 20/9*q + 8/81*q^2 - 1/9*s^2*q + 8/243*s^4:
gea4_f := (x, u) ->
  gea4_f0(gea4_s(x), gea4_q(u)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(gea4_f, rs, z, xs0, xs1, u0, u1):
