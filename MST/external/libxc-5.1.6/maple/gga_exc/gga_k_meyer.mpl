(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

(* eqn (5) *)
meyer_feta := y -> 1/2*(1 + (1 - y^2)*log((1 + y)/abs(1 - y))/(2*y)):

(* eqn (7) *)
meyer_lambda := y -> (1 - meyer_feta(y)) / (3 * y^2 * meyer_feta(y)):

(* enhancement factor from eqn (1) *)
meyer_f := x -> 1 + meyer_lambda(X2S*x/6)*x^2/(8*K_FACTOR_C):

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(meyer_f, rs, zeta, xs0, xs1):
