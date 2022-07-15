(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_x_rpw86_params
$include "gga_x_pw86.mpl"

lv_alpha := 0.02178:
lv_beta  := 1.15:
lv_muLV   := 0.8491/9:

lv_f0 := s ->
   + (1 + lv_muLV*s^2)/(1 + lv_alpha*s^6)
   + lv_alpha*s^6*pw86_f0(s)/(lv_beta + lv_alpha*s^6):

lv_f  := x -> lv_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(lv_f, rs, z, xs0, xs1):
