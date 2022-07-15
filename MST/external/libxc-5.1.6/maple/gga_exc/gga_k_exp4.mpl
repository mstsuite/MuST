(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

exp4_a1 := 199.81:
exp4_a2 := 4.3476:
exp4_c1 := 0.8524:
exp4_c2 := 1.2264:

# This is Eq. (40) of the paper.
exp4_f0 := s -> exp4_c1*(1 - exp(-exp4_a1*s^2)) + exp4_c2*(1 - exp(-exp4_a2*s^4)):
exp4_f  := x -> exp4_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(exp4_f, rs, zeta, xs0, xs1):
