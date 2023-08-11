(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

cac  := 1.467:
mtau := 4.5:
bcgp_pt := t -> t*sqrt((mtau + t)/(mtau + cac*t)):
arg := 0.5:
crg := 0.16667:
drg := 0.29633:
fbeta_modRasoltGeldart := rs -> (1.0 + arg*rs*(1+crg*rs)) / (1.0 + arg*rs*(1+drg*rs)):

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

(* override definition of tp *)
tp := (rs, z, xt) -> bcgp_pt(tt(rs, z, xt)):

(* override definition of mbeta *)
mbeta := (rs, t) -> params_a_beta*fbeta_modRasoltGeldart(rs):
