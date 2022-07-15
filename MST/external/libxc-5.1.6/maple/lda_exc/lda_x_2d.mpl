(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$define xc_dimensions_2d

ax := -4/3*sqrt(2)/Pi:

f := (rs, z) -> ax*f_zeta_2d(z)/rs:
