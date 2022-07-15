(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

(* in the paper we have beta_a = 0.066725 *)
beta_a := 0.066724550603149220:
beta_b := 0.1:
beta_c := 0.1778:

(* we redefine beta here *)
(* this is the Hu and Langreth expression *)
mbeta := (rs, t) -> beta_a*(1 + beta_b*rs)/(1 + beta_c*rs):
