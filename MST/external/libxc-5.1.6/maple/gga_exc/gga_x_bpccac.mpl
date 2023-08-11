(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_x_pbe_tca_params
$include "gga_x_pbe.mpl"

$define gga_x_pw91_params
$include "gga_x_pw91.mpl"

bpccac_malpha :=  1:
bpccac_mbeta  := 19:

bpccac_fab := x -> 1/(1 + exp(-bpccac_malpha*(x - bpccac_mbeta))):
bpccac_f   := x -> (1 - bpccac_fab(x))*pbe_f(x) + bpccac_fab(x)*pw91_f(x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(bpccac_f, rs, z, xs0, xs1):
