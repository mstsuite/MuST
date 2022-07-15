(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

params_a_alpha := 1:
$include "lda_x.mpl"

beta := rs -> (9*Pi/4)^(1/3)/(rs*M_C):
phi  := rs -> 1 - 1.5*(sqrt(1 + beta(rs)^2)/beta(rs) - arcsinh(beta(rs))/beta(rs)^2)^2:

f    := (rs, z) -> f_lda_x(rs, z)*phi(rs):
