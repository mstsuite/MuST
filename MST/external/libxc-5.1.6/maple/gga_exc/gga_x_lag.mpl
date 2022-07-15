(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

lag_a1 :=   0.041106:
lag_a2 :=   2.626712:
lag_a3 :=   0.092070:
lag_a4 :=   0.657946:

lag_f0 := s-> lag_a1 * s^lag_a2/(1 + lag_a3 * s^lag_a2)^lag_a4:
lag_f  := x-> lag_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(lag_f, rs, zeta, xs0, xs1):
