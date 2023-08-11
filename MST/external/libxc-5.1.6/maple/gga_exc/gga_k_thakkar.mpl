(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

thakkar_f0 := x -> 1 + 0.0055*x^2/(1 + 0.0253*x*arcsinh(x)):
thakkar_f1 := x -> -0.072*x/(1 + 2*4^(1/3)*x):

thakkar_f := x -> thakkar_f0(x) + thakkar_f1(x):

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(thakkar_f, rs, zeta, xs0, xs1):
