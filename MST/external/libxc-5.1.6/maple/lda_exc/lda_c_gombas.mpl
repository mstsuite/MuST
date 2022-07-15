(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

(* eq 26 *)
a1 := -0.0357:
a2 :=  0.0562:
b1 := -0.0311:
b2 :=  2.39:

(* eq 25 *)
f := (rs, zeta) -> a1/(1 + a2*rs/RS_FACTOR)
  + b1*log((rs/RS_FACTOR + b2)/(rs/RS_FACTOR)):
