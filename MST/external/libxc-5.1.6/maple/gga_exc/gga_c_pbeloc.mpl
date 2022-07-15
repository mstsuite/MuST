(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

pbeloc_b0 := 0.0375:
pbeloc_a  := 0.08:

(* we redefine beta here *)
mbeta := (rs, t) -> pbeloc_b0 + pbeloc_a*t^2*(1 - exp(-rs^2)):
