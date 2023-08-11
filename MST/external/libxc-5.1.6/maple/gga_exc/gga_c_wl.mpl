(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

a := -0.74860:
b :=  0.06001:
c :=  3.60073:
d :=  0.90000:

f_num := (z, xt) -> sqrt(1 - z^2)*(a + b*xt):
f_den := (rs, xs0, xs1) -> c + d*(xs0 + xs1) + rs:

f := (rs, zeta, xt, xs0, xs1) -> f_num(zeta, xt)/f_den(rs, xs0, xs1):
