(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_x_pbe_sol_params
$include "gga_x_pbe.mpl"

q2d_cc := 100:
q2d_c1 := 0.5217:

q2d_f1 := s -> pbe_f0(s)*(q2d_cc - s^4) + q2d_c1*s^3.5*(1 + s^2):
q2d_f2 := s -> q2d_cc + s^6:
q2d_f  := x -> q2d_f1(X2S*x)/q2d_f2(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(q2d_f, rs, zeta, xs0, xs1):
