(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

a :=  0.0311:
b := -0.048:
c :=  0.009:
d := -0.017:

f := (rs, zeta) -> a*log(rs) + b + c*rs*log(rs) + d*rs:
