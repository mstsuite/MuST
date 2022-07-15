(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_x_pbe_sol_params
$include "gga_x_pbe.mpl"

q1d_a := 0.06525:

q1d_f1 := s -> pbe_f0(s) + (s^2 + s^4)/(1 + s^4 + s^6)*(-pbe_f0(s)*s^2 + q1d_a):
q1d_f  := x -> q1d_f1(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(q1d_f, rs, zeta, xs0, xs1):
