(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

rge2_kappa := 0.8040:

rge2_den := s -> rge2_kappa + 1*MU_GE*s^2 + MU_GE^2*s^4/rge2_kappa:
rge2_f0  := s -> 1 + rge2_kappa * (1 - rge2_kappa/rge2_den(s)):
rge2_f   := x -> rge2_f0(X2S * x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(rge2_f, rs, zeta, xs0, xs1):
