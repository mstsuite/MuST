(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_ssb_sw_params *params;

  assert(p->params != NULL);
  params = (gga_x_ssb_sw_params * )(p->params);
*)

ssb_sw_f0 := s -> params_a_A
   + params_a_B*s^2/(1 + params_a_C*s^2)
   - params_a_D*s^2/(1 + params_a_E*s^4):
ssb_sw_f := x -> ssb_sw_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(ssb_sw_f, rs, zeta, xs0, xs1):
