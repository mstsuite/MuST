(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_tau_hcth_params *params;

  assert(p->params != NULL);
  params = (mgga_x_tau_hcth_params * ) (p->params);
*)

hcth_coeff_a := [0, 1, 0, -2, 0, 1]:

(* Equation (29) *)
hcth_gamX := 0.004:
hcth_ux   := x -> hcth_gamX*x^2/(1 + hcth_gamX*x^2):

hcth_gxl  := x -> add(params_a_cx_local [i]*hcth_ux(x)^(i-1), i=1..4):
hcth_gxnl := x -> add(params_a_cx_nlocal[i]*hcth_ux(x)^(i-1), i=1..4):

hcth_f    := (x, u, t) -> hcth_gxl(x) + hcth_gxnl(x)*mgga_series_w(hcth_coeff_a, 6, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(hcth_f, rs, z, xs0, xs1, u0, u1, t0, t1):