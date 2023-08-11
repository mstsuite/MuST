(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_x_tm.mpl"

reg_c := 3:
reg_d := 1.475:

(* text after equation 11 *)
regtm_f := (a, p) -> (1 - a)^3/(1 + (reg_d*a)^2)^(3/2) * exp(-reg_c*p):

(* equation 11 *)
regtm_zp := (a, p) -> 1 / (1 + 3/5*(a/(p+regtm_f(a,p)))):

(* equation 13 *)
regtm_w := zp -> (zp^2 + 3*zp^3) / (1 + zp^3)^2:

(* Collect all the pieces together *)
tm_w := (x, t) -> regtm_w(regtm_zp(tm_alpha(x, t), tm_p(x))):
