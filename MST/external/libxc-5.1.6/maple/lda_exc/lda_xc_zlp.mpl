(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

a0 := 0.93222*RS_FACTOR:
kk := 9.47362e-3*RS_FACTOR:

f := (rs, zeta) -> -a0*(1 - kk*log(1 + rs/kk)/rs)/rs:
