(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

lp90_c0 := 0.80569:
lp90_d0 := 3.0124e-3:
lp90_k  := 4.0743e-3:

(* Equation (60) *)
lp90_f := (rs, z, xt, us0, us1) ->
  - (lp90_c0 + lp90_d0*t_vw(z, xt, us0, us1))/(rs/RS_FACTOR + lp90_k):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  lp90_f(rs, z, xt, us0, us1):
