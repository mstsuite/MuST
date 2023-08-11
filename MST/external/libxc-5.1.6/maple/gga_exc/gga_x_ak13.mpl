(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_ak13_params *params;

  assert(p->params != NULL);
  params = (gga_x_ak13_params * )(p->params);
*)

ak13_f0 := s -> 1 + params_a_B1*s*log(1 + s) + params_a_B2*s*log(1 + log(1 + s)):
ak13_f  := x -> ak13_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(ak13_f, rs, zeta, xs0, xs1):
